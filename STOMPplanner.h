#ifndef STOMPPLANNER_H
#define STOMPPLANNER_H

#include "../../../libLAML/include/Distribution.h"
#include "../../../libLAML/include/Matlab.h"
#include "../../../libLAML/include/DenseMatrix.h"
#include "sdfenv.h"
#include <Eigen/Dense>
#include <vector>
#include "visualPath3D.h"
#include <iostream>
#include "hapticXremoteVel.h"
#include "stompInterrupt.h"
#include "spline.h" //cubic插值
#include <numeric>
#include <algorithm>
#include <nav_msgs/Odometry.h>
#include "Bspline.h"
#include <optional>
#include "processElevationMap.h"
#include <set>
#include <string>

using namespace std;

class STOMPplanner_base{

public:
    static const string obstacleConstriant;
    static const string footholdNumConstriant;
    static const string shortPathConstriant;

    STOMPplanner_base
    (int nSamples,int kPath,int dim,const vector<Eigen::MatrixXd>& DiscretePoint_loc
        ,const sdfenv & Env,const grid_map::GridMap & footholdEnv,const vector<string>& vecConstriant)
    :nSamples(nSamples),kPath(kPath),R_libLAML(nSamples-2,nSamples-2),Rinv_libLAML(nSamples-2,nSamples-2)
    ,M_libLAML(nSamples-2,nSamples-2),dim(dim)
    ,Env(Env),DiscretePoint_loc(DiscretePoint_loc)
    ,footholdEnv(footholdEnv)
    {
        for (auto & ele:vecConstriant){
            Constriants_.insert(ele);
        }
        Precompute();
    }

    const sdfenv & readsdf()const{
        return Env;
    }

    const grid_map::GridMap & readFootholdEnv()const{
        return footholdEnv;
    }

    int readnSamples()const{
        return nSamples;
    }

    int readDim()const{
        return dim;
    }

    virtual void visualThetaSeq(const Eigen::MatrixXd & thetaSeq,Eigen::Vector3d colorRGB)const = 0;// 可视化角度序列，可能是直接可视化为一条路径，也可能可视化为一组六足机器人，由派生类决定

    Eigen::MatrixXd getSolverInitThetaSeq(){ //得到solver的初始的角度序列
        if (!initthetaSeq_.has_value()){
            cout<<"错误调用"<<endl;
            throw exception();
        }
        return initthetaSeq_.value();
    }

    Eigen::MatrixXd solver(Interruptbase& interrupt,const int iteUpper=20){
        for(auto & ele:DiscretePoint_loc){
            if (!Env.checkValidLoc(TransLoc2XYZ(ele)[0])){
                cout<<"输入的待规划初始序列超出了范围"<<endl;
                throw exception();
                // return Eigen::MatrixXd{};
            }
        }

        
        // if (DiscretePoint_loc.size()!=2){ //计算示教或者rrt的路径的path
        //     Eigen::MatrixXd thetaseqBeginning(readDim(),DiscretePoint_loc.size());
        //     for(auto & ele:DiscretePoint_loc){
        //         thetaseqBeginning InverseCompute(ele)
        //     }
        // }

        Eigen::MatrixXd initThetaSeq = compute_initThetaSeq(DiscretePoint_loc); //初始得到的是角度序列
        Eigen::MatrixXd thetaSeqCur = initThetaSeq;
        initthetaSeq_ = thetaSeqCur;
        vector<Eigen::MatrixXd> TrajN(kPath,Eigen::MatrixXd(dim,nSamples));
        vector<Eigen::MatrixXd> ek(dim,Eigen::MatrixXd(kPath,nSamples)); //ek[0]是kPath行nSamples列的double
        // cout<<thetaSeqCur<<endl;
        // cout<<R<<endl;
        // cout<<"111111111211"<<endl;
        // cout<<thetaSeqCur.row(0)<<endl;
        // cout<<"222141"<<endl;
        // cout<<thetaSeqCur.row(1)<<endl;
        double Qtheta = TrajsumCost(thetaSeqCur); 
        cout<<Qtheta<<endl;
        // cout<<"11111111111"<<endl;
        double QthetaOld = 0;
        Eigen::MatrixXd pathCost(nSamples,kPath);

        int ite = 0;
        while(abs(Qtheta-QthetaOld)>0.0001){
            if (interrupt.isInterrupt()){ //由人以某种方式打断
                break;
            }
            ite++;
            if (ite>iteUpper){
                break;
            }

            QthetaOld = Qtheta;
            cout<<ite<<" "<<Qtheta<<endl;

            visualThetaSeq(thetaSeqCur,Eigen::Vector3d(0,0,1));
            
            stompCompute_NoisyTraj(thetaSeqCur,TrajN,ek);

            for(int i = 0;i<TrajN.size();i++){
                pathCost.col(i) = stompCompute_Cost(TrajN[i]);
            }
            Eigen::MatrixXd pathE = stompCompute_ELambda(pathCost);
            Eigen::MatrixXd pathProb = pathE; //nsamples行 kpath列
            for(int i = 0;i<pathProb.rows();i++){
                double Sum = pathProb.row(i).sum();
                for(int j = 0;j<pathProb.cols();j++){
                    if (Sum == 0){
                        pathProb(i,j) = 0;
                    }
                    else{
                        pathProb(i,j) /= Sum; 
                    }
                }
            }
            
            for(int dim_i = 0;dim_i<dim;dim_i++){
                Eigen::MatrixXd epsilon = ek[dim_i];
                Eigen::VectorXd dtheta(nSamples-2);
                Eigen::MatrixXd tmp = pathProb.transpose(); //tmp是kpath行，nsamples列
                for(int i = 0;i<tmp.rows();i++){
                    for(int j = 0;j<tmp.cols();j++){
                        tmp(i,j) *= epsilon(i,j); 
                    }
                }
                for(int i = 1;i<dtheta.size()-1;i++){
                    dtheta(i) = tmp.col(i).sum();
                }
                if (pathCost.sum()==4){
                    dtheta = Eigen::VectorXd::Zero(nSamples-2);
                }
                dtheta = M * dtheta;
                thetaSeqCur.block(dim_i,1,1,nSamples-2) += dtheta.transpose();
            }
            Qtheta = TrajsumCost(thetaSeqCur); 
        }
        cout<<"STOMP solver finish"<<endl;
        // visualThetaSeq(thetaSeqCur,Eigen::Vector3d(1,0,0),"STOMPcurpath");
        return thetaSeqCur;
    }

    const vector<Eigen::MatrixXd>& read_DiscretePoint_loc()const{
        return DiscretePoint_loc;
    }
protected:
    set<string> Constriants_;
private:
    const int nSamples;
    Eigen::MatrixXd R;
    Eigen::MatrixXd Rinv;
    Eigen::MatrixXd M;
    const int kPath; //每次迭代生成的随机路径
    vector<Eigen::MatrixXd> DiscretePoint_loc; //输入的初始的路径离散点 是有限个机器人的整体配置

    std::optional<Eigen::MatrixXd> initthetaSeq_{};
    const int dim; //有多少个变量待规划

    const sdfenv & Env; //sdf图 这是障碍物的环境，大片的完全不可落足的障碍物区域
    const grid_map::GridMap & footholdEnv; //这是可落足的环境，包括几何的离散可落足点和带有物理属性的

    double safeHold = 5; //表明安全裕度 越大越保守

    DenseMatrix R_libLAML;
    DenseMatrix Rinv_libLAML;
    DenseMatrix M_libLAML; //为了能够用到这个包的其他函数

    //输入几个离散的轨迹点，然后看是用cubic插值还是用直线直接相连
    virtual Eigen::MatrixXd  compute_initThetaSeq(vector<Eigen::MatrixXd> DiscretePoint_loc) = 0;

    //theta是dim维向量，表示规划的各个关节的取值 需要重写这个 映射到机器人的整体配置，像2D路径就是映射到机体位置,还有就是各个位置的防止碰撞的包络半径
    //就比如说对六足机器人整体来说，theta就是18维向量；假如把六足整体看成24个球包络而成的，一条腿由4个球包络，那么就有24个坐标点和24个半径，分别存储在loc和radius里
    virtual void ForwardCompute(const Eigen::VectorXd & theta,Eigen::MatrixXd & loc,Eigen::VectorXd & radius)const = 0;

    //从机器人的整体配置映射到关节序列
    virtual Eigen::VectorXd InverseCompute(const Eigen::MatrixXd  & point)const = 0;

    //知道机器人上一个配置和当前配置，得到配置各处的速度值 一个配置占据一个vectorXd的位置
    virtual Eigen::VectorXd velocityCompute(const Eigen::MatrixXd & locLast,const Eigen::MatrixXd & locNow)const = 0;

    //知道机器人的整体配置，得到各配置在世界坐标系下的坐标(x,y,z)或(x,y)，以便在sdf图中进行索引 这也相当于是一个解析loc的协议
    virtual vector<Eigen::VectorXd> TransLoc2XYZ(const Eigen::MatrixXd & loc)const = 0;

    //由障碍物产生的cost 仅针对路径上的某一点来说的 不是针对路径序列来说的
    double stompCost_obstacle_singlePoint(const Eigen::MatrixXd & loc,const Eigen::VectorXd & radius,const Eigen::VectorXd &velocity)const{
        Eigen::VectorXd costTmp = radius;
        vector<Eigen::VectorXd> vecLocXYZ = TransLoc2XYZ(loc);
        double cost = 0;
        for(int i = 0;i<costTmp.size();i++){
            if (!Env.getsdfValue(vecLocXYZ[i]).has_value()){
                costTmp(i) = 0;
            }
            else{   
                costTmp(i) += (safeHold-Env.getsdfValue(vecLocXYZ[i]).value());
                costTmp(i) = max(costTmp(i),0);
                costTmp(i) *= velocity(i);
            }
            
            cost+=costTmp(i);
        }
        return cost;
    }

    //仅计算由障碍产生的cost,是针对路径序列来计算的
    Eigen::VectorXd stompCost_obstacle(const Eigen::MatrixXd& thetaSeq)const{
        Eigen::MatrixXd loc;
        Eigen::VectorXd radius;
        Eigen::VectorXd velocity;
        Eigen::MatrixXd locLast{};
        Eigen::VectorXd result(nSamples);
        for (int j = 0;j<thetaSeq.cols();j++){
            ForwardCompute(thetaSeq.col(j),loc,radius);
            // cout<<loc<<" "<<radius<<endl;

            if (locLast.size()==0){
                Eigen::VectorXd tmpzero(radius.size());
                tmpzero = Eigen::VectorXd::Zero(radius.size());
                velocity = tmpzero; //第一次velocity全0
            }
            else{
                velocity = velocityCompute(locLast,loc);
            }
            
            locLast = loc;

            result(j) = stompCost_obstacle_singlePoint(loc,radius,velocity);
        }
        return result;
    }

    //还可以加其他的cost 比如对于六足来说落足点什么的
    virtual Eigen::VectorXd other_Cost(const Eigen::MatrixXd & thetaSeq)const = 0;

    //计算thetaSeq表示的路径的分值 包括障碍物及其他约束 返回值是nSamples维的向量
    Eigen::VectorXd stompCompute_Cost(const Eigen::MatrixXd& thetaSeq)const{
        Eigen::VectorXd result(thetaSeq.cols());
        for(int i = 0;i<result.size();i++){
            result(i) = 0;
        }
        if (Constriants_.find(STOMPConstriantName::obstacleConstriant)!=Constriants_.end()){
            result+=stompCost_obstacle(thetaSeq);
        }
        result+=other_Cost(thetaSeq);
        return result;
        // return stompCost_obstacle(thetaSeq) + other_Cost(thetaSeq);
    }

    void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)const {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols() - 1;
        
        if( colToRemove < numCols ) {
            matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
            matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
        }
        
        matrix.conservativeResize(numRows,numCols);
    }


    //计算某条轨迹上的所有点的代价和 还需要带上光滑项
    double TrajsumCost(const Eigen::MatrixXd& thetaSeq)const{
        Eigen::VectorXd Costi = stompCompute_Cost(thetaSeq);

        // cout<<Costi<<endl;

        Eigen::MatrixXd thetaSeqSub = thetaSeq;
        RemoveColumn(thetaSeqSub,0);
        RemoveColumn(thetaSeqSub,nSamples-1);
        // cout<<thetaSeqSub<<endl;

        double result = 0;

        // for(int i = 0;i<Costi.size();i++){
        //     std::cout<<Costi(i)<<" ";
        // }
        // std::cout<<std::endl;
        result+=Costi.sum()+  0.5*(thetaSeqSub*R*thetaSeqSub.transpose()).sum();
        return result;
    }


    //以thetaSeq为基准，生成kPath条随机路径 thetaSeq是dim行nSamples列
    //ek.size = dim ek[0]是kPath行nSamples列的double
    //Traj_N.size = kpath Traj_N[0]是dim行nSamples列
    void stompCompute_NoisyTraj(const Eigen::MatrixXd& thetaSeq,vector<Eigen::MatrixXd>&Traj_N,vector<Eigen::MatrixXd> & ek) {
        if (ek.size()!=dim || Traj_N.size()!=kPath){
            throw exception();
            std::cout<<"ek.size()!=dim || Traj_N.size()!=kPath"<<std::endl;
            // ek = vector<Eigen::MatrixXd>(dim);
            // Traj_N = vector<Eigen::MatrixXd>(kPath);
        }
    
        DenseMatrix mu(1,Rinv.cols());
        for(int i = 0;i<mu.getColumnDimension();i++){
            mu.setEntry(0,i,0);
        }
        Eigen::MatrixXd eigenTmp(kPath,nSamples);
        eigenTmp = Eigen::MatrixXd::Zero(kPath,nSamples);
        for(int i = 0;i<ek.size();i++){
            Matrix &tmp = mvnrnd(mu,Rinv_libLAML,kPath,"whatever");
            Eigen::MatrixXd tmpblock(kPath,nSamples-2); ;
            TranslibLAML2Eigen(tmp,tmpblock);
            eigenTmp.block(0,1,kPath,nSamples-2) = tmpblock;
            ek[i] = eigenTmp;
            // cout<<"1111111111111111111111111111111111-----------------------------------------------------"<<endl;
            // cout<<ek[i]<<endl;
        }

        for (int i = 0;i<kPath;i++){
            for (int dim_i = 0;dim_i<dim;dim_i++){
                Traj_N[i].row(dim_i)= thetaSeq.row(dim_i) + ek[dim_i].row(i);
            }
        }
    }

    //exp(xxx) pathCost是nsamples行kPath列
    Eigen::MatrixXd stompCompute_ELambda(const Eigen::MatrixXd & pathCost)const {
        double h = 10;  //是一个常值，这个常值在论文中被设置为了10
        Eigen::VectorXd maxS(nSamples);
        Eigen::VectorXd minS(nSamples);
        for(int i = 0;i<nSamples;i++){
            maxS(i) = pathCost.row(i).maxCoeff();
            minS(i) = pathCost.row(i).minCoeff();
        }
        Eigen::MatrixXd tmp = pathCost;
        for (int j = 0;j<tmp.cols();j++){
            tmp.col(j) = tmp.col(j)-minS;
        }
        Eigen::VectorXd minusMaxMin = maxS-minS;
        for(int j=0;j<tmp.cols();j++){
            for(int i = 0;i<tmp.rows();i++){
                if (minusMaxMin(i)!=0){
                    tmp(i,j) = exp(-h*tmp(i,j)/minusMaxMin(i));
                }
                else{
                    tmp(i,j) = 0;
                }
            }
        }
        return tmp;
    }

    void TransEigen2libLAML(const Eigen::MatrixXd & mat,DenseMatrix & mat_libLAML)const{
        for(int i = 0;i<mat.rows();i++){
            for(int j = 0;j<mat.cols();j++){
                mat_libLAML.setEntry(i,j,mat(i,j));
            }
        }
    }
    void TranslibLAML2Eigen(Matrix & mat_libLAML,Eigen::MatrixXd & mat)const{
        for(int i = 0;i<mat_libLAML.getRowDimension();i++){
            for(int j = 0;j<mat_libLAML.getColumnDimension();j++){
                mat(i,j) = mat_libLAML.getEntry(i,j);
            }
        }
    }

    void Precompute(){
        Eigen::MatrixXd A(nSamples,nSamples-2);
        A = Eigen::MatrixXd::Zero(nSamples,nSamples-2); //有这个初始化很重要，不然会出现莫名其妙的错误
        for(int j=0;j<A.cols();j++){
            A(j+0,j) = 1;
            A(j+1,j) = -2;
            A(j+2,j) = 1;
        }
        R =  A.transpose()*A;
        // cout<<A<<endl;
        // cout<<A.transpose()<<endl;
        Rinv = R.inverse();

        M = 1.0/nSamples * Rinv;
        Eigen::MatrixXd maxRinv(1,Rinv.cols());
        for(int j=0;j<Rinv.cols();j++){
            maxRinv(0,j) = Rinv.col(j).maxCoeff();
        }
        for (int i = 0;i<M.rows();i++){
            for(int j = 0;j<M.cols();j++){
                M(i,j) = M(i,j)/maxRinv(0,j);
            }
        }
        Rinv = Rinv/Rinv.sum();

        TransEigen2libLAML(R,R_libLAML);
        TransEigen2libLAML(Rinv,Rinv_libLAML);
        TransEigen2libLAML(M,M_libLAML);
    }
};

const string STOMPplanner_base::obstacleConstriant = STOMPConstriantName::obstacleConstriant;
const string STOMPplanner_base::footholdNumConstriant = STOMPConstriantName::footholdNumConstriant;
const string STOMPplanner_base::shortPathConstriant = STOMPConstriantName::shortPathConstriant;

class STOMPHexapodPath:public STOMPplanner_base{
public:
    STOMPHexapodPath
                        (int nSamples,int kPath,int dim,const sdfenv & Env
                        ,const grid_map::GridMap & footholdEnv
                        ,const vector<Eigen::MatrixXd>& DiscretePoint_loc
                        ,double PlannerHeiht = 0.5,double collisionRadius = 20,double scaleLower=0,double scaleUpper=0.3
                        ,vector<string> vecConstriant = vector<string>({STOMPConstriantName::obstacleConstriant}) )
        :scaleLower(scaleLower),scaleUpper(scaleUpper)
        ,STOMPplanner_base(nSamples,kPath,dim,DiscretePoint_loc,Env,footholdEnv,vecConstriant)
        ,collisionRadius(collisionRadius)
        ,PlannerHeiht(PlannerHeiht)
        {
            
        }

    //用于将角度序列转换到待跟踪的路径序列,对于六足的路径规划来说，由于需要兼容旧代码的缘故，需要做这个事情
    vector<geometry_msgs::Point> thetaSeq2PointSeq(const Eigen::MatrixXd & thetaSeq){
        Eigen::MatrixXd loc;
        Eigen::VectorXd radius;
        vector<Eigen::VectorXd> vecPoint;
        for(int j = 0;j<thetaSeq.cols();j++){
            ForwardCompute(thetaSeq.col(j),loc,radius);
            vecPoint.push_back(TransLoc2XYZ(loc)[0]);
        }
        vector<geometry_msgs::Point> result;
        for(auto &ele:vecPoint){
            geometry_msgs::Point pntTmp;
            pntTmp.x = ele(0);
            pntTmp.y = ele(1);
            pntTmp.z = ele(2);
            result.push_back(pntTmp);
        }
        return result;
    }

private:
    // Eigen::Vector2d startPoint; //世界坐标系下的X Y
    // Eigen::Vector2d endPoint;

    // Eigen::Vector2d startTheta; //对应的经过缩放的值
    // Eigen::Vector2d endTheta; 

    double scaleLower;
    double scaleUpper; 

    double collisionRadius; //能够完全包络机体的圆的半径大小

    double PlannerHeiht; //定义STOMP2D规划器的规划平面高度 默认为0.5m，但是考虑到感知到的实际障碍物在elevationmapping中大多数比0.5小，因此会有一个用户自定义的值

    template<typename Type>
    vector<Type> SampleN(vector<Type> vecPoint,int num){
        if (num>=vecPoint.size()){
            return vecPoint;
        }
        int sampleinterval = num; //采样sampleinterval个点
        int numindex = 0;
        for(int i = 0;i<vecPoint.size();i+=vecPoint.size()/sampleinterval){
            swap(vecPoint[numindex],vecPoint[i]);
            numindex++;
        }
        if (numindex<vecPoint.size()){
            swap(vecPoint[numindex],vecPoint[vecPoint.size()-1]);
            numindex++;
        }
        vecPoint.erase(vecPoint.begin()+numindex,vecPoint.end());
        return vecPoint;
    }

    virtual Eigen::MatrixXd compute_initThetaSeq(vector<Eigen::MatrixXd> DiscretePoint_loc)override{
        int type = 1; //如果不是直线那么就默认采用B样条插值
        Eigen::MatrixXd initThetaSeq(readDim(),readnSamples()); 
        if (DiscretePoint_loc.size()==2){ //这种情况直接按直线插值处理
            Eigen::VectorXd startTheta = InverseCompute(DiscretePoint_loc[0]);
            Eigen::VectorXd endTheta = InverseCompute(DiscretePoint_loc[1]);
            Eigen::VectorXd step = 1.0/(readnSamples()-1) * (endTheta-startTheta);
            for(int j = 0;j<readnSamples();j++){
                initThetaSeq.col(j) = startTheta + j*step;
            }
        }
        else if (type == 1){
            
            vector<geometry_msgs::Point> vecPoint;
            for(int i = 0;i<DiscretePoint_loc.size();i++){
                geometry_msgs::Point pnt;
                pnt.x = DiscretePoint_loc[i](0);
                pnt.y = DiscretePoint_loc[i](1);
                pnt.z = 0;
                vecPoint.push_back(pnt);
            }
            
            int sampleTimes = 6;
            vecPoint = SampleN(vecPoint,sampleTimes);

            Bspline bspline(vecPoint);
            vector<geometry_msgs::Point> output = bspline.solve();
            output = SampleN(output,readnSamples());
            while(output.size()!=readnSamples()){
                output.pop_back();
            }
            if (output.size()!=readnSamples()){
                cout<<"output.size()!=readnSamples()"<<output.size()<<endl;
                throw exception();
            }
            // cout<<"123124125135126"<<endl;
            DiscretePoint_loc.clear();
            for(int i = 0;i<output.size();i++){
                // cout<<Eigen::Vector3d(output[i].x,output[i].y,0)<<endl;
                DiscretePoint_loc.push_back(Eigen::Vector3d(output[i].x,output[i].y,0));
            }
            cout<<DiscretePoint_loc.size()<<endl;
            for(int i = 0;i<DiscretePoint_loc.size();i++){
                initThetaSeq.col(i) = InverseCompute(DiscretePoint_loc[i]);
            }

        }
        else {
            
            //认为一般示教的路径朝向y轴方向
            auto prev = accumulate(DiscretePoint_loc.begin(),DiscretePoint_loc.begin()+DiscretePoint_loc.size()/2,Eigen::Vector3d(0,0,0));
            prev /= DiscretePoint_loc.size()/2;
            auto after = accumulate(DiscretePoint_loc.begin()+DiscretePoint_loc.size()/2,DiscretePoint_loc.end(),Eigen::Vector3d(0,0,0));
            after /= DiscretePoint_loc.size() - DiscretePoint_loc.size()/2;
            int sigh ;
            if (after(1)>prev(1)){
                sigh = 1;
            }
            else{
                sigh = -1;
            }
            
            vector<double> X,Y;
            double cur;
            if (sigh==1){
                cur = numeric_limits<double>::lowest();
            }
            else{
                cur = numeric_limits<double>::max();
            }
            for(int i = 0;i<DiscretePoint_loc.size();i++){
                if (sigh==1 && DiscretePoint_loc[i](1)>cur){
                    cur = DiscretePoint_loc[i](1);
                    X.push_back(sigh*DiscretePoint_loc[i](1)); //因为必须是严格增的
                    Y.push_back(DiscretePoint_loc[i](0));
                }
                else if (sigh==-1 && DiscretePoint_loc[i](1)<cur){
                    cur = DiscretePoint_loc[i](1);
                    X.push_back(sigh*DiscretePoint_loc[i](1)); //因为必须是严格增的
                    Y.push_back(DiscretePoint_loc[i](0));
                }
                // else {
                //     cout<<"erro"<<endl;
                //     throw exception();
                // }
                
            }


            int sampleinterval = 5; //采样sampleinterval个点，然后进行cubic插值
            X = SampleN(X,sampleinterval);
            Y = SampleN(Y,sampleinterval);

            tk::spline s(X,Y);
            DiscretePoint_loc.clear();
            for(double i = X[0];i<*(X.end()-1);i+=(*(X.end()-1)-X[0])/readnSamples()){
                DiscretePoint_loc.push_back(Eigen::Vector3d(s(i),sigh*i,0));
                // cout<<"x spline:"<<i<<" "<<s(i)<<endl;
            }

            if (DiscretePoint_loc.size()>readnSamples()){ //输入的DiscretePoint_loc一定是不重复的点集序列
                while (DiscretePoint_loc.size()!=readnSamples()){ //去除多余的
                    DiscretePoint_loc.pop_back();
                }
            }
            cout<<DiscretePoint_loc.size()<<endl;
            
            
            for(int i = 0;i<DiscretePoint_loc.size();i++){
                initThetaSeq.col(i) = InverseCompute(DiscretePoint_loc[i]);
            }
        }

        return initThetaSeq; // Eigen::MatrixXd initThetaSeq; //初始的角度序列，dim行nSamples列
    }

    //输入是一个列向量point,表示在世界坐标系下的位置
    virtual Eigen::VectorXd InverseCompute(const Eigen::MatrixXd& point)const override {
        Eigen::Vector2d Point2D {point(0),point(1)};
        Eigen::VectorXd result(2);
        result = (scaleUpper-scaleLower)*(Point2D-readsdf().getLowerLeft());
        result(0) /=  readsdf().getSize()(0);
        result(1) /=  readsdf().getSize()(1);
        return result;
    }


    virtual void ForwardCompute(const Eigen::VectorXd & theta,Eigen::MatrixXd & loc,Eigen::VectorXd & radius)const override{
        Eigen::MatrixXd result(2,1);
        result = Eigen::Vector2d{theta(0)-scaleLower,theta(1)-scaleLower};
        result = 1.0/(scaleUpper-scaleLower)*result;
        result(0) *= readsdf().getSize()(0);
        result(1) *= readsdf().getSize()(1);
        result += readsdf().getLowerLeft();
        loc = result;
        Eigen::VectorXd radiusret(1);
        radiusret(0) = collisionRadius;
        radius = radiusret;
    }

    virtual Eigen::VectorXd velocityCompute(const Eigen::MatrixXd & locLast,const Eigen::MatrixXd & locNow)const override{
        Eigen::VectorXd result(1); //对于六足机器人的路径规划来说，一个值就够了
        Eigen::Vector2d vel = locNow-locLast;
        vel(0) *=vel(0);
        vel(1) *=vel(1);
        result(0) = sqrt(vel.sum());
        return result;
    }

    virtual vector<Eigen::VectorXd> TransLoc2XYZ(const Eigen::MatrixXd & loc)const override{
        return vector<Eigen::VectorXd>{Eigen::Vector3d{loc(0,0),loc(1,0),PlannerHeiht}};
    }

    Eigen::VectorXd computeFootholdNumvalue(const Eigen::MatrixXd & thetaseq)const{ //计算沿路径的落足点密度，以便衡量路径的好坏
        const grid_map::GridMap & footholdmap = readFootholdEnv();

        double detectSize = 2; //探测周围两米的范围内的落足点数目 为了编程方便，范围是一个正方形范围
        
        Eigen::VectorXd result(thetaseq.cols());
        result = Eigen::VectorXd::Zero(thetaseq.cols());
        for(int j = 0;j<thetaseq.cols();j++){
            Eigen::VectorXd PointLocNow = getPosfromTheta(thetaseq.col(j));

            grid_map::Index luIndex;
            grid_map::Index rdIndex;
            footholdmap.getIndex(grid_map::Position(PointLocNow(0)-detectSize/2,PointLocNow(1)-detectSize/2),luIndex);
            footholdmap.getIndex(grid_map::Position(PointLocNow(0)+detectSize/2,PointLocNow(1)+detectSize/2),rdIndex);
            for(int i_ = min(luIndex(0),rdIndex(0));i_<max(luIndex(0),rdIndex(0));i_++){
                for(int j_ = min(luIndex(1),rdIndex(1));j_<max(luIndex(1),rdIndex(1));j_++){
                    grid_map::Index indextmp(i_,j_);
                    
                    if (!std::isnan(footholdmap.at(layerName::geometryFootholdLayer,indextmp))){
                        result(j)++; //只要在给定的范围内有一个几何可落足点，就进行计数
                    }
                }
            }
        }

        double MAXfootholdNum = readnSamples()* pow(detectSize/footholdmap.getResolution(),2);
        for(int i = 0;i<result.size();i++){
            result(i) = MAXfootholdNum-result(i);
            result(i) *= 0.0001;
        }
        return result;
    }

    Eigen::VectorXd computeFootholdPhysicalValue(const Eigen::MatrixXd & thetaseq)const{ //计算沿路径的落足点密度，以便衡量路径的好坏
        const grid_map::GridMap & footholdmap = readFootholdEnv();

        double detectSize = 2; //探测周围两米的范围内的落足点 为了编程方便，范围是一个正方形范围
        
        Eigen::VectorXd result(thetaseq.cols());
        result = Eigen::VectorXd::Zero(thetaseq.cols());
        for(int j = 0;j<thetaseq.cols();j++){
            Eigen::VectorXd PointLocNow = getPosfromTheta(thetaseq.col(j));

            grid_map::Index luIndex;
            grid_map::Index rdIndex;
            footholdmap.getIndex(grid_map::Position(PointLocNow(0)-detectSize/2,PointLocNow(1)-detectSize/2),luIndex);
            footholdmap.getIndex(grid_map::Position(PointLocNow(0)+detectSize/2,PointLocNow(1)+detectSize/2),rdIndex);
            for(int i_ = min(luIndex(0),rdIndex(0));i_<max(luIndex(0),rdIndex(0));i_++){
                for(int j_ = min(luIndex(1),rdIndex(1));j_<max(luIndex(1),rdIndex(1));j_++){
                    grid_map::Index indextmp(i_,j_);
                    
                    if (!std::isnan(footholdmap.at(layerName::physicalFootholdLayer,indextmp))){
                        result(j)+=footholdmap.at(layerName::physicalFootholdLayer,indextmp); //只要在给定的范围内有一个物理落足点，就加
                    }
                }
            }
        }
        //物理特征的最大值是1 归一化的 0-1的取值范围 1是最好的物理特征

        double MAXfootholdphysical = readnSamples()* pow(detectSize/footholdmap.getResolution(),2);
        for(int i = 0;i<result.size();i++){
            result(i) = MAXfootholdphysical-result(i);
            result(i) *= 0.0001;
        }
        return result;
    }

    Eigen::VectorXd getPosfromTheta(const Eigen::VectorXd & theta)const{ //从角度值得到位置
        Eigen::MatrixXd loc;
        Eigen::VectorXd radius;
        ForwardCompute(theta,loc,radius);
        return TransLoc2XYZ(loc)[0];
    }

    Eigen::VectorXd computeRouteLength(const Eigen::MatrixXd & thetaseq)const{ //计算路径长度 不能让路径太长了 还是有点作用的
        double ref = (getPosfromTheta(thetaseq.col(0))-getPosfromTheta(thetaseq.col(thetaseq.cols()-1))).norm();
        double sumPath = 0;
        for(int j = 0;j<thetaseq.cols()-1;j++){
            Eigen::VectorXd PointLocNow = getPosfromTheta(thetaseq.col(j));
            Eigen::VectorXd PointLocNext = getPosfromTheta(thetaseq.col(j+1));
            sumPath+=(PointLocNext-PointLocNow).norm();
        }
        sumPath = sumPath -1*ref;
        Eigen::VectorXd result(thetaseq.cols());
        result = Eigen::VectorXd::Zero(thetaseq.cols());
        for(int j=0;j<result.size();j++){
            result(j) = sumPath/readnSamples();
            result(j)*=0.1; //这个值可以决定有多大程度上听最短路径的话 如果这个值很大的话，基本上就是直线了，否则很小的话就基本上听其他的代价项了
        }
        return result;
    }

    virtual Eigen::VectorXd other_Cost(const Eigen::MatrixXd & thetaSeq)const override{ //暂时先设置成均为0
        // Eigen::VectorXd result = computeRouteLength(thetaSeq);
        //  Eigen::VectorXd result = computeFootholdNumvalue(thetaSeq);
        //  Eigen::VectorXd result = computeRouteLength(thetaSeq)+computeFootholdNumvalue(thetaSeq); //加上这个很容易就不收敛了

        Eigen::VectorXd result(thetaSeq.cols());
        for(int i = 0;i<result.size();i++){
            result(i) = 0;
        }
        if (Constriants_.find(STOMPConstriantName::footholdNumConstriant)!=Constriants_.end()){
            result+=computeFootholdNumvalue(thetaSeq);
        }
        if (Constriants_.find(STOMPConstriantName::shortPathConstriant)!=Constriants_.end()){
            result+=computeRouteLength(thetaSeq);
        }
        if (Constriants_.find(STOMPConstriantName::physicalConstriant)!=Constriants_.end()){
            result+=computeFootholdPhysicalValue(thetaSeq);
        }
        return result;
        
    }

    virtual void visualThetaSeq(const Eigen::MatrixXd & thetaSeq,Eigen::Vector3d colorRGB)const override{
        // cout<<topicname<<"sssssssssssssssssssssssssssss"<<endl;
        static visualPath3D visPath("world","STOMPcurpath"); //用来可视化的
        vector<Eigen::Vector3d> Pointsvis;
        Eigen::MatrixXd thetaseqvis = thetaSeq;

        //可视化的
        for(int j = 0;j<thetaseqvis.cols();j++){
            Eigen::VectorXd thetaNow = thetaseqvis.col(j);
            Eigen::MatrixXd tmploc ;
            Eigen::VectorXd tmp1;
            ForwardCompute(thetaNow,tmploc,tmp1);
            vector<Eigen::VectorXd> locXYZ = TransLoc2XYZ(tmploc);
            if(!readsdf().checkValidLoc(locXYZ[0])){
                cout<<locXYZ[0]<<endl;
                Eigen::VectorXd thetaNow = thetaseqvis.col(j);
                Eigen::MatrixXd tmploc ;
                Eigen::VectorXd tmp1;
                ForwardCompute(thetaNow,tmploc,tmp1);
                vector<Eigen::VectorXd> locXYZ = TransLoc2XYZ(tmploc);
                getchar(); //这个只是调试用的，没有什么作用 过一段时间就删了吧
            }
            Pointsvis.push_back(locXYZ[0]);
        }
        visPath.visual(Pointsvis,colorRGB);
    }
};



#endif