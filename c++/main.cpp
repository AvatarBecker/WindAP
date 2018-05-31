s#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include "opencv2/opencv.hpp"
#include "types.hpp"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>

#include <iomanip>

#include <chrono>

//#define USE_PID
#define USE_SCHEDULED_PID
//#define USE_THE_FORCE_LUKE
//#define USE_CAMACHO
//#define USE_ZHU_ZHI
//#define PLANT_IDENT
//#define TRACKBARS_ON

using namespace cv;
using namespace std;

const unsigned cCam = 1;
const double cBar_Dist_Real = .415;

chrono::milliseconds cT(70);  // dt in ms
double cT_sec = double(cT.count())/1000;

const int cPWM_Res = 1023; // pwm resolution
const double cMax_PWM = 1023;
const double cMin_PWM = 0;
const double cMin_Y = 0.05;
const double cMax_Y = 0.36;

#if defined USE_PID

    double cKP = 1120;
    double cKI = 600;
    double cKD = 1650;

    int cKP_temp = 0;   // to receive from trackbar...
    int cKI_temp = 0;
    int cKD_temp = 0;

    // convert to digital:
    double a = cKP + cKI*cT_sec/2 + cKD/cT_sec;
    double b = cKI*cT_sec/2-2*cKD/cT_sec - cKP;
    double c = cKD/cT_sec;

#elif defined USE_SCHEDULED_PID

    double cY_LOW  = 0.111;
    double cKP_LOW = 1097;
    double cKI_LOW = 874;
    double cKD_LOW = 1535;

    double cY_HIGH  = 0.263;
    double cKP_HIGH = 1300;
    double cKI_HIGH = 700;
    double cKD_HIGH = 962;

    double cKP(double y)
    {
        double m = (cKP_HIGH-cKP_LOW)/(cY_HIGH-cY_LOW);
        return (m*(y-cY_LOW)+cKP_LOW);
    }

    double cKI(double y)
    {
        double m = (cKI_HIGH-cKI_LOW)/(cY_HIGH-cY_LOW);
        return (m*(y-cY_LOW)+cKI_LOW);
    }

    double cKD(double y)
    {
        double m = (cKD_HIGH-cKD_LOW)/(cY_HIGH-cY_LOW);
        return (m*(y-cY_LOW)+cKD_LOW);
    }

    // convert to digital:
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

#elif defined USE_THE_FORCE_LUKE
    double cC2         = 221.65;            // fluid mechanics constant...
    double cTau        = .32;               // tau for the exp decay
    double cR          = exp(-cT_sec/cTau); // discrete exponential ratio
    double cVel_To_PWM = 583.24;            // Convert wind velocity to PWM (hoping to be linear..)
    double cG          = 9.78;

    int    cVel_To_PWM_temp    = 0;
#elif defined USE_CAMACHO
    double cA      = .125;
    double cLambda = .975;

    int cA_temp = 0;        // to receive from trackbar...
    int cLambda_temp = 0;
#elif defined USE_ZHU_ZHI
    double cTau    = 2;   // analog, first order decay constant (desired)
    double cLambda = .5;
#endif

const int cMin_Hue = 0;
const int cMax_Hue = 255;
const int cMin_xSaturation = 0;
const int cMax_Saturation = 255;
const int cMin_Value = 0;
const int cMax_Value = 255;

const int cBall_Min_Hue = 109;
const int cBall_Max_Hue = 133;
const int cBall_Min_Saturation = 100;
const int cBall_Max_Saturation = 190;
const int cBall_Min_Value = 70;
const int cBall_Max_Value = 160;

const int cBall_Min_Area = 900;

const int cBar_Min_Hue = 70;
const int cBar_Max_Hue = 100;
const int cBar_Min_Saturation = 100;
const int cBar_Max_Saturation = 255;
const int cBar_Min_Value = 48;
const int cBar_Max_Value = 142;

const Scalar cBall_Min(cBall_Min_Hue, cBall_Min_Saturation, cBall_Min_Value);
const Scalar cBall_Max(cBall_Max_Hue, cBall_Max_Saturation, cBall_Max_Value);

const Scalar cBar_Min(cBar_Min_Hue, cBar_Min_Saturation, cBar_Min_Value);
const Scalar cBar_Max(cBar_Max_Hue, cBar_Max_Saturation, cBar_Max_Value);

const double  cCenter_Neighborhood_Threshold = 70.0;
const double  cBar_Min_M00                   = 200.0;

const string cWindow_Original = "WindAPP";
const string cWindow_Binary = "Binary";

// how to return multiple objects? pass a vector of them and push back each one found?
//int DetectObject(Mat image, Scalar thres_min, Scalar thres_max);

bool IsAreaLess(Moments lhs, Moments rhs)
{
    return lhs.m00 < rhs.m00;
}

double Distance(Point P0, Point P1)
{
    return sqrt((double)pow(P0.x-P1.x,2) + (double)pow(P0.y-P1.y,2));
}

double sgn(double T)
{
    return T/abs(T);
}

int main()
{
    // CONTROL STUFF
    int   pwm = 923;
    double y_r = 0.25;
    int y_r_temp = 0.0;
    double y_r_temp_double = 0.0;

    #if defined PLANT_IDENT

        ofstream out_file("windap_data5.dat");
        int pwm_change_counter = 0;
        int max_pwm = 900;
        int min_pwm = 765;

    #elif defined USE_PID || defined USE_SCHEDULED_PID

        vector<double> e{0.0, 0.0, 0.0};
        vector<double> u{0.0, 0.0};

    #elif defined USE_CAMACHO
        // for recursive least squares (RLS)
        Mat M     = 100*(Mat::eye(4,4,CV_64F));
        Mat theta = (Mat_<double>(4,1) << 0.003, 0.006, 0.75, 0.05);
        Mat gain  = (Mat_<double>(4,1) << 0.0, 0.0, 0.0, 0.0);     // the gamma

        double rls_error = 0.0;     // difference between predict and actual y

//        Mat y = (Mat_<double>(2,1) << 0.0, 0.0);       // ditto, dy
        Mat y = (Mat::zeros(3,1,CV_64F));  // ditto, dy
        //Mat du = (Mat_<double>(4,1) << 0.0, 0.0, 0.0);  // ditto, du
        Mat u = (Mat::zeros(7,1,CV_64F));  // ditto, du
        Mat e  = (Mat_<double>(3,1) << 0.0, 0.0, 0.0);       // ditto, u

        Mat x  = (Mat_<double>(4,1) << 0.0, 0.0, 0.0, 0.0);  // memory vector for one k

        double b0 = 0.0;
        double b1 = 0.0;
        double a1 = 0.0;
        double a2 = 0.0;

        double g0 = 0.0;
        double g2 = 0.0;
        double g3 = 0.0;

    #elif defined USE_ZHU_ZHI

        unsigned cLength_Theta = 3;

        // for recursive least squares (RLS)
        Mat M     = 10000*(Mat::eye(cLength_Theta, cLength_Theta, CV_64F));
        Mat theta = 10000*(Mat::ones(cLength_Theta, 1, CV_64F));
        Mat gain  = (Mat::zeros(cLength_Theta, 1, CV_64F));     // the gamma

        double rls_error = 500;     // difference between predict and actual u

        Mat u = 1000*(Mat::ones(cLength_Theta, 1, CV_64F));
        Mat y = .03*(Mat::ones(cLength_Theta, 1, CV_64F));
        Mat e = -.15*(Mat::ones(cLength_Theta, 1, CV_64F));

        Mat x = .03*(Mat::ones(cLength_Theta, 1, CV_64F));  // memory vector for one k

        Mat y_bar = (Mat::zeros(cLength_Theta, 1, CV_64F));
        Mat u_bar = (Mat::zeros(cLength_Theta+1, 1, CV_64F));

        double d0 = 1000;
        double d1 = 1000;
        double d2 = 1000;
        double a  = 2/(cTau*cT_sec);

        double KP = 0;
        double KI = 0;
        double KD = 0;

    #endif

    // SOCKET STUFF
    struct sockaddr_in myaddr;
    int sock;
    string message = "100";
    memset(&myaddr, 0, sizeof(myaddr));
    myaddr.sin_family = AF_INET;
    myaddr.sin_addr.s_addr = htonl(INADDR_ANY);
    myaddr.sin_port = htons(0);
    inet_pton(AF_INET, "192.168.4.1", &myaddr.sin_addr.s_addr);
//    inet_pton(AF_INET, "127.0.0.1", &myaddr.sin_addr.s_addr); //send to local address (receive in with terminal for debugging)
    myaddr.sin_port = htons(3031);

    if( (sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0 )
    {
        perror("Failed to create socket");
    }

    // OPENCV STUFF
    Mat frame_raw;
    Mat frame_flipped;
    Mat frame_hsv;
    Mat frame_bin_ball;
    Mat frame_bin_bar;
    //Mat frame_edge_ball;

    vector<Moments>::iterator largest_contour;

    VideoCapture camera;
    camera.open(cCam);
    namedWindow(cWindow_Original, CV_WINDOW_AUTOSIZE);
    //namedWindow(cWindow_Binary, CV_WINDOW_AUTOSIZE);

    #if defined TRACKBARS_ON
        #if defined USE_CAMACHO
            cvCreateTrackbar("A", "WindAPP", &cA_temp, 10000);
            cvCreateTrackbar("lambda", "WindAPP", &cLambda_temp, 10000);
        #elif defined USE_PID
            cvCreateTrackbar("KP", "WindAPP", &cKP_temp, 10000);
            cvCreateTrackbar("KI", "WindAPP", &cKI_temp, 10000);
            cvCreateTrackbar("KD", "WindAPP", &cKD_temp, 10000);
        #elif defined USE_THE_FORCE_LUKE
            cvCreateTrackbar("KD", "WindAPP", &cVel_To_PWM_temp, 10000);
        #endif
    #endif

    cvCreateTrackbar("setpoint", "WindAPP", &y_r_temp, 10000);



    while(1)
    {
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        chrono::steady_clock::time_point finish;
        chrono::milliseconds elapsed_time;
        chrono::milliseconds remaining_time;

        double   bar_dist_px = 20;  // distance between bars as represented in the image (in pixels)
        double   factor_px2real = 1.0;
        double   ball_height_px = 0;
        double   ball_height_real = 0;

        vector<Moments> moments_ball;
        vector<Moments> moments_bar;
        Point   center_ball(0,0);
        vector<Point> centers_bars;
        vector< vector<Point> > contours;
        vector< Vec4i > hierarchy;

        camera >> frame_raw;
        //flip(frame_raw, frame_raw, 1);

        cvtColor(frame_raw, frame_hsv, COLOR_BGR2HSV, 0);

        // Find bars and print it to the raw frame:
        inRange(frame_hsv, cBar_Min, cBar_Max, frame_bin_bar);

        erode(frame_bin_bar, frame_bin_bar, getStructuringElement(MORPH_RECT, Size(3,3)));
        dilate(frame_bin_bar, frame_bin_bar, getStructuringElement(MORPH_RECT, Size(7,7)));

        findContours(frame_bin_bar, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);

        if( hierarchy.size() > 0 )
        {
            // for each contour, find moments, compare area.
            for(int i=0; i<contours.size(); ++i)
            {
                Moments moment_bar_temp;

                drawContours(frame_raw, contours, i, Scalar(0, 255, 255));
                moment_bar_temp = moments(contours[i],true);

                if( moment_bar_temp.m00 > cBar_Min_M00 )
                {
                    moments_bar.push_back(moment_bar_temp);
                }
            }

            Point center_temp_prev(10000, 10000);   //
            for(int i = 0; i < moments_bar.size(); ++i)
            {
                Point center_temp;

                center_temp.x = moments_bar[i].m10/moments_bar[i].m00;
                center_temp.y = moments_bar[i].m01/moments_bar[i].m00;

                if( Distance(center_temp, center_temp_prev) > cCenter_Neighborhood_Threshold )
                {
                    centers_bars.push_back(center_temp);
                    //cout << "centers " << i << ": " << center_temp.x << ", " << center_temp.y << endl;
                    circle(frame_raw, center_temp, 2, Scalar(255, 0, 255), 2);

                    center_temp_prev = center_temp;
                }
            }
        }

        // Find ball and print it to the raw frame:
        inRange(frame_hsv, cBall_Min, cBall_Max, frame_bin_ball);

        erode(frame_bin_ball, frame_bin_ball, getStructuringElement(MORPH_ELLIPSE, Size(3,3)));
        dilate(frame_bin_ball, frame_bin_ball, getStructuringElement(MORPH_ELLIPSE, Size(7,7)));

        findContours(frame_bin_ball, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);

        // If ball contours were found...
        if( hierarchy.size() > 0 )
        {
            // for each contour, find moments, compare area.
            for(int i=0; i<contours.size(); i++)
            {
                moments_ball.push_back(moments((Mat)contours[i]));
            }

            largest_contour = max_element(moments_ball.begin(), moments_ball.end(), IsAreaLess);
            int index_largest_contour = distance(moments_ball.begin(), largest_contour) - 1;
            index_largest_contour = (index_largest_contour >= 0) ? index_largest_contour : 0;

            center_ball.x = largest_contour->m10/largest_contour->m00;
            center_ball.y = largest_contour->m01/largest_contour->m00;

            // drawObject(x, y, frame_raw);
            circle(frame_raw, center_ball, 2, Scalar(255, 0, 255), 2);
            //drawContours(frame_raw, contours, index_largest_contour, Scalar(0, 255, 0));
        }

        // Iff two bars were found, than we can trust the measure and use it for controlling the plant
        if( centers_bars.size() == 2 )
        {
            bar_dist_px = Distance(centers_bars[0], centers_bars[1]);
            factor_px2real = cBar_Dist_Real/bar_dist_px;

            ball_height_px = Distance(center_ball, centers_bars[0]);
            ball_height_real = factor_px2real*ball_height_px;

            y_r_temp_double = double(y_r_temp);
            y_r = .03 + .34*y_r_temp_double/10000;  // update setpoint

            stringstream height_str;
            height_str << "h  = " << fixed << setprecision(1) << 100*ball_height_real;
            putText(frame_raw, height_str.str(), Point(10, 20),  FONT_HERSHEY_SIMPLEX, .5, Scalar(0, 255, 255), 1);

            stringstream sp_str;
            sp_str << "sp = " << fixed << setprecision(1) << 100*y_r;
            putText(frame_raw, sp_str.str(), Point(10, 45),  FONT_HERSHEY_SIMPLEX, .5, Scalar(0, 255, 255), 1);

            cout << pwm << " " << 100*ball_height_real << " " << 100*y_r << endl;

            #if defined USE_PID

                e[2] = e[1];
                e[1] = e[0];
                e[0] = y_r - ball_height_real;

                u[1] = u[0];
                u[0] = u[1] + a*e[0] + b*e[1] + c*e[2];

                pwm = int(u[0]);

                cout << "pwm pid: " << pwm << endl;

                if( u[0] > cMax_PWM ) pwm = cMax_PWM;
                if( u[0] < cMin_PWM ) pwm = cMin_PWM;

                #if defined TRACKBARS_ON
                    double cKP_double = double(cKP_temp);
                    double cKI_double = double(cKI_temp);
                    double cKD_double = double(cKD_temp);

                    cKP = 3000*double(cKP_double)/10000;
                    cKI = 3000*double(cKI_double)/10000;
                    cKD = 3000*double(cKD_double)/10000;

                    cout << "Kp: " << cKP << "; Ki: " << cKI << "; Kd: " << cKD << "\t e: " << e[0] << endl << endl;
                #endif

                // convert to digital:
                a = cKP + cKI*cT_sec/2 + cKD/cT_sec;
                b = -cKP + cKI*cT_sec/2 - 2*cKD/cT_sec;
                c = cKD/cT_sec;

            #elif defined USE_SCHEDULED_PID

                double y = ball_height_real;

                // convert to digital:
                a = cKP(y_r) + cKI(y_r)*cT_sec/2 + cKD(y_r)/cT_sec;
                b = -cKP(y_r) + cKI(y_r)*cT_sec/2 - 2*cKD(y_r)/cT_sec;
                c = cKD(y_r)/cT_sec;

                e[2] = e[1];
                e[1] = e[0];
                e[0] = y_r - y;

                u[1] = u[0];
                u[0] = u[1] + a*e[0] + b*e[1] + c*e[2];

                pwm = int(u[0]);

                cout << "pwm pid: " << pwm << endl;

                if( u[0] > cMax_PWM ) pwm = cMax_PWM;
                if( u[0] < cMin_PWM ) pwm = cMin_PWM;

                cout << "Kp: " << cKP(y_r) << "; Ki: " << cKI(y_r) << "; Kd: " << cKD(y_r) << "\t e: " << e[0] << endl << endl;

            #elif defined USE_THE_FORCE_LUKE

                if( ball_height_real > cMax_Y )
                {
                    pwm = 300;
                    cout << " max!" << endl << endl;
                }
                else if( ball_height_real < cMin_Y )
                {
                    pwm = 1000;
                    cout << " min!" << endl << endl;
                }
                else
                {
                    double e = y_r - ball_height_real;
                    double arg_temp = 1 - cG/cC2 - e*cR/(cC2*cTau*cTau) + e*cR*cR/(cTau*cTau);
                    pwm = cVel_To_PWM*(1+sgn(arg_temp)*sqrt(abs(arg_temp)));
                    cout << "e: " << e << "; pwm cont: " << pwm << endl << endl;
                }

                if( pwm > cMax_PWM ) pwm = cMax_PWM;
                if( pwm < cMin_PWM ) pwm = cMin_PWM;

                #if defined TRACKBARS_ON
                    double cVel_To_PWM_double = double(cVel_To_PWM_temp);

                    cVel_To_PWM = 500 + 1000*cVel_To_PWM_double/10000;

                    cout << endl << "cVel_To_PWM: " << endl << cVel_To_PWM << endl << endl;
                #endif

            #elif defined USE_CAMACHO

                y.at<double>(2) = y.at<double>(1);
                y.at<double>(1) = y.at<double>(0);
                y.at<double>(0) = ball_height_real;

//                if( y.at<double>(0) > .04 && y.at<double>(0) < .36 )
//                {
                    // --- RLS ---
                    // mount memory vector
                    x.at<double>(0) = u.at<double>(u.rows-2);
                    x.at<double>(1) = u.at<double>(u.rows-1);
                    x.at<double>(2) = -y.at<double>(y.rows-2);
                    x.at<double>(3) = -y.at<double>(y.rows-1);

                    Mat_<double> temp(1,1);

                    temp = ( x.t()*theta );
                    rls_error = y.at<double>(0) - temp.at<double>(0);

                    temp = ( x.t()*M*x );
                    gain  = M*x/(cLambda + temp.at<double>(0));

                    theta = theta + gain*rls_error;
                    M = (M - gain*x.t()*M)/cLambda;
                    // --- \RLS ---
//                }

                e.at<double>(2) = e.at<double>(1);
                e.at<double>(1) = e.at<double>(0);
                e.at<double>(0) = y_r - y.at<double>(0);

                // calculation of the control law
                b0 = theta.at<double>(0);
                b1 = theta.at<double>(1);
                a1 = theta.at<double>(2);
                a2 = theta.at<double>(3);

                g0 = pow(cA,u.rows-1)*(1-cA)/(b0*cA+b1);
                double arg_temp = a1*a1-4*a2;
                g2 = 0.5*( a1 + sgn(arg_temp)*sqrt( abs(arg_temp) ) );
                g3 = a2/g2;

                // store previous
                for(int i = 0; i < u.rows-1; ++i )
                {
                    u.at<double>(i+1) = u.at<double>(i);
                }

//                if( y.at<double>(0) > .36 )
//                {
//                    u.at<double>(0) = cMin_PWM;
//                }
//                else if( y.at<double>(0) < .04 )
//                {
//                    u.at<double>(0) = cMax_PWM;
//                }
//                else
//                {
                    u.at<double>(0) = 1000000*g0*( e.at<double>(0) - (g2+g3)*e.at<double>(1) + g2*g3*e.at<double>(2) );  // set control variable
//                }

                cout << "pwm cont: " << u.at<double>(0) << endl;

                // limiter
                if( u.at<double>(0) > cMax_PWM ) u.at<double>(0) = cMax_PWM;
                if( u.at<double>(0) < cMin_PWM ) u.at<double>(0) = cMin_PWM;

                pwm = (int)u.at<double>(0);

                #if defined TRACKBARS_ON
                    double cA_double = double(cA_temp);
                    double cLambda_double = double(cLambda_temp);

                    cA     = .05+2*double(cA_double)/10000;
                    cLambda = .1 + .9*double(cLambda_double)/10000;
                #endif

                cout << endl << "theta: " << endl << theta << endl;
                cout << endl << "a: " << cA << ", lambda: " << cLambda << endl << endl;

            #elif defined USE_ZHU_ZHI

                y.at<double>(2) = y.at<double>(1);
                y.at<double>(1) = y.at<double>(0);
                y.at<double>(0) = ball_height_real;

                y_bar.at<double>(2) = y_bar.at<double>(1);
                y_bar.at<double>(1) = y_bar.at<double>(0);
                y_bar.at<double>(0) = a*( y.at<double>(1) - y.at<double>(0) );

                u_bar.at<double>(3) = u_bar.at<double>(2);
                u_bar.at<double>(2) = u_bar.at<double>(1);
                u_bar.at<double>(1) = u_bar.at<double>(0);
                u_bar.at<double>(0) = u.at<double>(1) - u.at<double>(3);

//                if( y.at<double>(0) > .04 && y.at<double>(0) < .36 )
//                {
                    // --- RLS ---
                    // mount memory vector
                    x.at<double>(0) = y_bar.at<double>(0);
                    x.at<double>(1) = y_bar.at<double>(1);
                    x.at<double>(2) = y_bar.at<double>(2);
                    //cout << "x: " << endl << x << endl;

                    Mat_<double> temp(1,1);

                    temp = ( x.t()*theta );
                    //cout << "temp: " << endl << temp << endl;
                    rls_error = u_bar.at<double>(0) - temp.at<double>(0);
                    //cout << "rls_error: " << endl << rls_error << endl;

                    temp = ( x.t()*M*x );
                    //cout << "temp: " << endl << temp << endl;
                    gain  = (M*x)/( cLambda + temp.at<double>(0) );
                    //cout << "gain: " << endl << gain << endl;

                    theta = theta + gain*rls_error;
                    //cout << "theta: " << endl << theta << endl;
                    M = (M - gain*x.t()*M)/cLambda;
                    //cout << "M: " << endl << M << endl;
                    // --- \RLS ---
//                }

                e.at<double>(2) = e.at<double>(1);
                e.at<double>(1) = e.at<double>(0);
                e.at<double>(0) = y_r - y.at<double>(0);

                // calculation of the control law
                d0 = theta.at<double>(0);
                d1 = theta.at<double>(1);
                d2 = theta.at<double>(2);

                KP = -( d1 + 2*d2 );
                KI = -cT_sec*( d1 + 2*d2 )/( d0 + d1 + d2 );
                KD = -cT_sec*d1/( d1 + 2*d2 );

                // store previous
                u.at<double>(2) = u.at<double>(1);
                u.at<double>(1) = u.at<double>(0);

//                if( y.at<double>(0) > .36 )
//                {
//                    u.at<double>(0) = cMin_PWM;
//                }
//                else if( y.at<double>(0) < .04 )
//                {
//                    u.at<double>(0) = cMax_PWM;
//                }
//                else
//                {
                    u.at<double>(0) = u.at<double>(1) + d0*e.at<double>(0) + d1*e.at<double>(1) + d2*e.at<double>(2);  // set control variable
//                }

                cout << "pwm cont: " << u.at<double>(0) << "; e: " << e.at<double>(0) << endl;

                // limiter
                if( u.at<double>(0) > cMax_PWM ) u.at<double>(0) = cMax_PWM;
                if( u.at<double>(0) < cMin_PWM ) u.at<double>(0) = cMin_PWM;

                pwm = (int)u.at<double>(0);

                cout << "Kp: " << KP << "; Ki: " << KI << "; Kd: " << KD << endl;
                cout << "d0: " << d0 << "; d1: " << d1 << "; d2: " << d2 << endl << endl;

            #elif defined PLANT_IDENT

                out_file << pwm << " " << ball_height_real << endl;

                if(pwm_change_counter > 17)
                {
//                    pwm = rand() % (cPWM_Res - cMin_PWM + 1) + cMin_PWM;

                    pwm = (pwm == min_pwm) ? max_pwm : min_pwm -= 2, min_pwm;
                    pwm_change_counter = 0;
                }
                else
                {
                    ++pwm_change_counter;
                }
            #endif

            // send PWM value
            int pwm_inv = cPWM_Res - pwm;   // inverted logic...
            string message = to_string(pwm_inv);
            sendto(sock, message.c_str(), message.size(), 0, (struct sockaddr *)&myaddr, sizeof(myaddr) );
        }
        else
        {
            cout << "Too many or too few reference bars detected!" << endl;
            cout << "Adjust lighting/filters." << endl << endl;

            putText(frame_raw, "h = ???", Point(10, 20),  FONT_HERSHEY_SIMPLEX, .5, Scalar(0, 255, 255), 1);
        }





        imshow(cWindow_Original, frame_raw);
        //imshow(cWindow_Binary, frame_bin_bar + frame_bin_ball);

        // Check loop time and try and keep it constant:
        finish = chrono::steady_clock::now();
//        elapsed_time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
        elapsed_time = chrono::duration_cast<chrono::milliseconds>(finish - start);

        remaining_time = cT - elapsed_time; // time left?

        cout << "dt = " << elapsed_time.count() << " ms; remaining = " << remaining_time.count() << endl;

        if( remaining_time.count() > 0 )
        {
            if(waitKey(remaining_time.count()) == 27) // wait for remaining_time. If ESC is pressed, break
            {
                break;
            }
        }
        else
        {
            cout << endl << "Loop is taking longer than " << cT.count() <<"ms! Close other programs or use GPU!" << endl;
        }
    }

    return 0;
}
