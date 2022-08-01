 <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <vector>
using namespace std;

//在rviz中画点 线
class visual_rviz{
public:
    visual_rviz(ros::NodeHandle& n,string topicName,string frame)
    :n_(n)
    ,topicName_(topicName)
    ,frame_(frame)
    ,marker_pub_(n.advertise<visualization_msgs::Marker>(topicName, 10))
    {
        //初始化
        points_.header.frame_id = line_strip_.header.frame_id = line_list_.header.frame_id = frame;
        points_.header.stamp = line_strip_.header.stamp = line_list_.header.stamp = ros::Time::now();
        points_.ns = line_strip_.ns = line_list_.ns = "points_and_lines";
        points_.action = line_strip_.action = line_list_.action = visualization_msgs::Marker::ADD;
        points_.pose.orientation.w = line_strip_.pose.orientation.w = line_list_.pose.orientation.w = 1.0;

        //分配3个id
        points_.id = 0;
        line_strip_.id = 1;
        line_list_.id = 2;

        //初始化形状
        points_.type = visualization_msgs::Marker::POINTS;
        line_strip_.type = visualization_msgs::Marker::LINE_STRIP;
        line_list_.type = visualization_msgs::Marker::LINE_LIST;
    }

    void addPoints(const vector<geometry_msgs::Point>& vec,Eigen::Vector2d scale,Eigen::Vector3d color){
        pointsStorage_.insert(pointsStorage_.end(),vec.begin(),vec.end());
        scalePoints_ = scale;
        vector<Eigen::Vector3d> vecColor(vec.size(),color);
        colorsPoints_.insert(colorsPoints_.end(),vecColor.begin(),vecColor.end());
    }
    void addPoints(const vector<vector<geometry_msgs::Point>>& vecVec,Eigen::Vector2d scale,vector<Eigen::Vector3d> colors){
        for (int i = 0;i<vecVec.size();i++){
            addPoints(vecVec[i],scale,colors[i%colors.size()]);
        }
    }
    void addPoints(const geometry_msgs::Point & point,Eigen::Vector2d scale,Eigen::Vector3d color){
        pointsStorage_.push_back(point);
        scalePoints_ = scale;
        colorsPoints_.push_back(color);
    }

    void addLineLists(const vector<geometry_msgs::Point>& vec,double scale,Eigen::Vector3d color){
        lineListsStorage_.insert(lineListsStorage_.end(),vec.begin(),vec.end());
        scaleLineLists_ = scale;
        vector<Eigen::Vector3d> vecColor(vec.size(),color);
        colorsLineLists_.insert(colorsLineLists_.end(),vecColor.begin(),vecColor.end());
    }
    void addLineLists(const vector<vector<geometry_msgs::Point>>& vecVec,double scale,vector<Eigen::Vector3d> colors){
        for (int i = 0;i<vecVec.size();i++){
            addLineLists(vecVec[i],scale,colors[i%colors.size()]);
        }
    }

    void visPoints(){
        visPoints(pointsStorage_,scalePoints_(0),scalePoints_(1),colorsPoints_);
        pointsStorage_.clear();
        colorsPoints_.clear();
    }
    void visLineLists(){
        visLineList(lineListsStorage_,scaleLineLists_,colorsLineLists_);
        lineListsStorage_.clear();
        colorsLineLists_.clear();
    }
    

private:
    ros::NodeHandle n_;
    string topicName_;
    string frame_;
    ros::Publisher marker_pub_;

    visualization_msgs::Marker points_, line_strip_, line_list_;

    vector<geometry_msgs::Point> pointsStorage_;
    vector<Eigen::Vector3d> colorsPoints_;
    Eigen::Vector2d scalePoints_; 

    vector<geometry_msgs::Point> lineListsStorage_;
    vector<Eigen::Vector3d> colorsLineLists_;
    double scaleLineLists_;



    void visPoints(const vector<geometry_msgs::Point>& vec,double scalex,double scaley,double r,double g,double b){
        points_.color.a = 1.0;
        points_.color.r = r;
        points_.color.g = g;
        points_.color.b = b;
        points_.points = vec;
        points_.scale.x = scalex;
        points_.scale.y = scaley;
        marker_pub_.publish(points_);
    }
    //把待可视化的展开成一串进行可视化
    void visPoints(const vector<geometry_msgs::Point>& vec,double scalex,double scaley,const vector<Eigen::Vector3d>&colors){
        std_msgs::ColorRGBA rgb;
        vector<std_msgs::ColorRGBA> vecrgb;
        for(auto & ele:colors){
            rgb.a = 1;
            rgb.r = ele(0);
            rgb.g = ele(1);
            rgb.b = ele(2);
            vecrgb.push_back(rgb);
        }
        points_.colors = vecrgb;
        points_.points = vec;
        points_.scale.x = scalex;
        points_.scale.y = scaley;
        marker_pub_.publish(points_);
    }
    //分层次的进行可视化 这个更直观更好用一些
    //vecVec.size()==colors.size()
    void visPoints(const vector<vector<geometry_msgs::Point>>& vecVec,double scalex,double scaley,const vector<Eigen::Vector3d>&colors){
        if (vecVec.size()!=colors.size()){
            cout<<"not vecVec.size()==colors.size()"<<endl;
            throw exception();
        }

        //展开成一维，然后进行可视化
        vector<Eigen::Vector3d> colors1dim;
        vector<geometry_msgs::Point> vecpoints1dim;
        for (int i = 0;i<vecVec.size();i++){
            vector<Eigen::Vector3d> colorTmp(vecVec[i].size(),colors[i]);
            colors1dim.insert(colors1dim.end(),colorTmp.begin(),colorTmp.end());
            vecpoints1dim.insert(vecpoints1dim.end(),vecVec[i].begin(),vecVec[i].end());
        }
        visPoints(vecpoints1dim,scalex,scaley,colors1dim);
    }

    void visLineList(const vector<geometry_msgs::Point>& vec,double scalex,double r,double g,double b){
        line_list_.color.a = 1.0;
        line_list_.color.r = r;
        line_list_.color.g = g;
        line_list_.color.b = b;
        line_list_.points = vec;
        line_list_.scale.x = scalex;
        marker_pub_.publish(line_list_);
    }
    void visLineList(const vector<geometry_msgs::Point>& vec,double scalex,const vector<Eigen::Vector3d>&colors){
        std_msgs::ColorRGBA rgb;
        vector<std_msgs::ColorRGBA> vecrgb;
        for(auto & ele:colors){
            rgb.a = 1;
            rgb.r = ele(0);
            rgb.g = ele(1);
            rgb.b = ele(2);
            vecrgb.push_back(rgb);
        }
        line_list_.colors = vecrgb;
        line_list_.points = vec;
        line_list_.scale.x = scalex;
        marker_pub_.publish(line_list_);
    }
    void visLineList(const vector<vector<geometry_msgs::Point>>& vecVec,double scalex,const vector<Eigen::Vector3d>&colors){
        if (vecVec.size()!=colors.size()){
            cout<<"not vecVec.size()==colors.size()"<<endl;
            throw exception();
        }

        //展开成一维，然后进行可视化
        vector<Eigen::Vector3d> colors1dim;
        vector<geometry_msgs::Point> vecpoints1dim;
        for (int i = 0;i<vecVec.size();i++){
            vector<Eigen::Vector3d> colorTmp(vecVec[i].size(),colors[i]);
            colors1dim.insert(colors1dim.end(),colorTmp.begin(),colorTmp.end());
            vecpoints1dim.insert(vecpoints1dim.end(),vecVec[i].begin(),vecVec[i].end());
        }
        visLineList(vecpoints1dim,scalex,colors1dim);
    }

    void visLineStrip(const vector<geometry_msgs::Point>& vec,double scalex,double r,double g,double b){
        line_strip_.color.a = 1.0;
        line_strip_.color.r = r;
        line_strip_.color.g = g;
        line_strip_.color.b = b;
        line_strip_.points = vec;
        line_strip_.scale.x = scalex;
        marker_pub_.publish(line_strip_);
    }
};
Footer
© 2022 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
