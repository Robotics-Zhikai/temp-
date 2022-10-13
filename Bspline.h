#ifndef BSPLINE_H
#define BSPLINE_H
#include "EISpider_description/planning.h"
#include"string.h"
#include<iostream>
#include<sstream>
#include<fstream>
#include <vector>
#include <geometry_msgs/Point.h>
using namespace std;

class Bspline
{
public:
    Bspline(const std::vector<geometry_msgs::Point> & inPoint):inPoint(inPoint){}
    std::vector<geometry_msgs::Point> solve(){
        return DrawBspline1(inPoint);
    }
private:
    const std::vector<geometry_msgs::Point> & inPoint;

    std::vector<geometry_msgs::Point> DrawBspline1(std::vector<geometry_msgs::Point> inPoint)
    {
        int num = inPoint.size();
        //系数矩阵对角列
        float *a = new float[num];
        //系数矩阵对角上列
        float *b = new float[num];
        //系数矩阵对角下列
        float *c = new float[num];

        //定义soluctionX、soluctionY为线性方程的解
        float *soluctionX = new float[num];
        float *soluctionY = new float[num];
        //定义dataX和dataY,用来存放inPoint里的X和Y坐标
        float *dataX = new float[num];
        float *dataY = new float[num];
        //定义controlPoint用来存放控制点
        geometry_msgs::Point *controlPoint = new geometry_msgs::Point[num + 4];
        //存放画线的两个使用点
        geometry_msgs::Point *lines = new geometry_msgs::Point[2];

        //初始化 a,b,c
        a[0] = 18;
        a[num - 1] = 18;

        for (int i = 1; i < num - 1; i++)
        {
            a[i] = 4;
        }

        for (int i = 1; i < num - 1; i++)
        {
            b[i] = 1;
            c[i] = 1;
        }

        c[num - 1] = -9;
        b[0] = -9;



        for (int i = 0; i < num; i++)
        {
            dataX[i] = 6.0f * inPoint[i].x;
            dataY[i] = 6.0f * inPoint[i].y;
        }

        dataX[0] *= 1.5f;
        dataY[0] *= 1.5f;
        dataX[num - 1] *= 1.5f;
        dataY[num - 1] *= 1.5f;

        //调用Matrix用追赶法求解线性方程
        Matrix(dataX, num, a, b, c, soluctionX);
        Matrix(dataY, num, a, b, c, soluctionY);

        controlPoint[num + 3].x = dataX[num - 1] / 9;
        controlPoint[num + 2].x = dataX[num - 1] / 9;
        controlPoint[0].x = dataX[0] / 9;
        controlPoint[1].x = dataX[0] / 9;


        for (int i = 0; i < num; i++)
        {
            controlPoint[i + 2].x = soluctionX[i];

        }

        controlPoint[num + 3].y = dataY[num - 1] / 9;
        controlPoint[num + 2].y = dataY[num - 1] / 9;
        controlPoint[0].y = dataY[0] / 9;
        controlPoint[1].y = dataY[0] / 9;


        for (int i = 0; i < num; i++)
        {
            controlPoint[i + 2].y = soluctionY[i];

        }
        std::vector<geometry_msgs::Point> oupPutPointList;

        // float interval = (num+1)/(float)outnum;
        
        for (int i = 0; i < num + 1; i++)
        {

            for (float u = 0.01; u <= 1; u += 0.01)
            {
                float b0 = 1.0f / 6 * (1 - u) * (1 - u) * (1 - u);
                float b1 = 1.0f / 6 * (3 * u * u * u - 6 * u * u + 4);
                float b2 = 1.0f / 6 * (-3 * u * u * u + 3 * u * u + 3 * u + 1);
                float b3 = 1.0f / 6 * u * u * u;
                geometry_msgs::Point dataP;
                dataP.x = (b0 * controlPoint[i].x + b1 * controlPoint[i + 1].x + b2 * controlPoint[i + 2].x + b3 * controlPoint[i + 3].x);
                dataP.y = (b0 * controlPoint[i].y + b1 * controlPoint[i + 1].y + b2 * controlPoint[i + 2].y + b3 * controlPoint[i + 3].y);
                oupPutPointList.push_back(dataP);
            }
        }
        delete[] a,b,c,soluctionX, soluctionY, dataX, dataY, controlPoint, lines;
        return oupPutPointList;
    }


    void Matrix(float * constantTerm, int num, float * m, float * n, float * k, float * solution)
    {
        //b为分解后的下三角矩阵的对角数组
        float *a = new float[num];
        //a为分解后的单位上三角矩阵的对角上方数组
        float *b = new float[num - 1];
        //c为分解后的单位上三角矩阵的对角上方数组
        float *c = new float[num];
        //x为求解过程中的间接解
        float *x = new float[num];
        int i;

        a[0] = m[0];
        b[0] = n[0] / a[0];

        //给分解后下三角矩阵的对角下方数组c赋值
        for (i = 1; i < num; i++)
        {
            c[i] = k[i];
        }


        //给分解后的单位上三角矩阵的对角上方数组a和分解后的单位上三角矩阵的对角上方数组b赋值
        for (i = 1; i < num - 1; i++)
        {
            a[i] = m[i] - k[i] * b[i - 1];
            b[i] = n[i] / a[i];

        }

        a[num - 1] = m[num - 1] - k[num - 1] * b[num - 2];
        //中间解x的初始值
        x[0] = constantTerm[0] / a[0];

        //给中间解赋值
        for (i = 1; i < num; i++)
        {
            x[i] = (constantTerm[i] - k[i] * x[i - 1]) / a[i];
        }

        //解出最终解
        solution[num - 1] = x[num - 1];

        for (i = num - 1; i > 0; i--)
        {
            solution[i - 1] = x[i - 1] - solution[i] * b[i - 1];
        }
        
        delete[] a,b,c,x;
    }

};






// int main( int argc, char** argv )
// {
//   ros::init(argc, argv, "threeGaitPlan");
//   ros::NodeHandle nh_;
//   ros::Duration(0.5).sleep();  
//   ros::Rate r(1000);

   
//     vector<vector<float>> obstacleList{
//         {7, 5, 1},
//         {5, 6, 2},
//         {5, 8, 2},
//         {5, 10, 2},
//         {9, 5, 2},
//         {11, 5, 2}
//     };

//     // ��ʼ���Ŀ���
//     node* startNode = new node(2.0, 2.0);
//     node* goalNode = new node(14.0, 9.0);

//     RRT rrt(startNode, goalNode, obstacleList, 0.5, 5);
//     vector<node*> aa = rrt.planning();
//     for(auto &ele:aa){
//       cout << ele->x << "," << ele->y << endl;
//     }
//     return 0;


// }

#endif