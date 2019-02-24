#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>//求逆等矩阵运算
#include "math.h"//指数运算
#include <random>
#include <vector>

using namespace Eigen;
using namespace std;

void norm_srand(double m,double sd,int num,vector<double> &N)//输入均值m，标准差sd
{   //采用Box-Muller得到服从正态分布的随机数
    double PI = 3.141592654;

    double U, V;
    int phase = 0;
    double z;
    srand((unsigned)time(NULL));
    for(int i=0;i<num;i++)
    {
        if(phase == 0)
        {
            U = rand() / (RAND_MAX + 1.0);
            V = rand() / (RAND_MAX + 1.0);
            z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
        }
        else
        {
            z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
        }
        phase = 1 - phase;
        N.push_back(z*sd+m);//将生成的均值1,标准差0的数转换成均值m,标准差sd
    }
}


int main()
{
    Vector4d x(0.0, 0.0, 0.0, 0.0);//状态变量  x位移，y位移，x速度，y速度
    Vector4d P0(1000.0, 1000.0, 1000.0, 1000.0);
    MatrixXd P(P0.asDiagonal());//由向量生成对角阵,行人的不确定性因素（先验估计协方差矩阵）
    double const dt = 0.1;//时间步长
    Matrix4d F;
    F << 1.0, 0.0, dt, 0.0,
            0.0, 1.0, 0.0, dt,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0;//定义状态转移矩阵
    Matrix<double, 2, 4> H;
    H << 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0;//定义观测矩阵

    //计算测量过程的协方差矩阵R和过程噪声的协方差矩阵Q
    double ra = pow(10.0, 2);
    Matrix<double, 2, 2> R;
    R << ra, 0.0,
            0.0, ra;//定义噪声的协方差矩阵
    double sv = 0.5;
    Vector4d G(0.5 * pow(dt, 2), 0.5 * pow(dt, 2), dt, dt);
    Matrix<double, 4, 4> Q;
    Q = G * G.transpose() * pow(sv, 2);

    //产生一些测量数据
    int num = 200;//测量数据的个数
    double vx = 20, vy = 10;

    vector<double> measurements_1,measurements_2,mx,my;
    norm_srand(0, 1, num, measurements_1);//正态分布的随机数vector
    norm_srand(0, 1, num, measurements_2);
    for(int j=0;j<measurements_1.size();j++)
    {
        mx.push_back(vx+measurements_1[j]);
        my.push_back(vy+measurements_2[j]);
        cout << mx[j] << " ";
    }
    cout << endl<<endl;

    //卡尔曼滤波-开始计算
    Matrix<double,2,2> S;
    Matrix<double,4,2> K;//卡尔曼增益
    Vector2d z;          //观测值
    cout<<"x位移\t"<<"y位移\t"<<"x速度\t"<<"y速度\t"<<endl;
    for(int i=0;i<mx.size();i++)
    {
        z << mx[i],my[i];                      //观测值
        x = x.transpose()*F;                   //根据物理方程的预测值
        P = F*P*F.transpose()+Q;               //预测值的误差
        S = H*P*H.transpose()+R;
        K = P*H.transpose()*S.inverse();       //计算卡尔曼增益

        x = x+K*(z-(H*x));                     //计算k+1时刻最优估计值

        P = ( MatrixXd::Identity(4,4)-(K*H))*P;//计算k+1时刻最优估计值的误差
        cout<<"第"<<i+1<<"次估计"<<endl;
        cout<<"速度观测值：\t"<<"        \t"<<z[0]<<"\t"<<z[1]<<endl;
        cout<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\t"<<x[3]<<endl<<endl;
    }
}