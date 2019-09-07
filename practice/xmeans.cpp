
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/bgsegm.hpp>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
#define Pi 3.1415
int K = 1;

Ptr<BackgroundSubtractor> pMOG; //MOG Background subtractor
//コンパイル: g++ TEST.cpp -I/usr/local/include/opencv2 -I/usr/local/include/opencv -L/usr/lib -lopencv_core -lopencv_imgcodecs -lopencv_highgui -lopencv_imgproc -lopencv_bgsegm

int xmeans(Mat &mat, Mat &orders, Mat &labels, int *K)
{
    double xu = 0.0;
    double yu = 0.0;
    double variancex = 0.0;
    double variancey = 0.0;
    double variancextmp = 0.0;
    double varianceytmp = 0.0;
    int Rn = mat.rows;
    printf("Rn = %d\n", Rn);
    int p = 2;
    int q = 2 * p;
    for (int i = 0; i < Rn; i++)
    {
        xu += mat.at<float>(i, 0);
        yu += mat.at<float>(i, 1);
        if (i % 1000 == 0)
        {
            //printf(" x = %f, y = %f\n", mat.at<float>(i, 0), mat.at<float>(i, 1));
        }
    }
    xu = xu / Rn;
    yu = yu / Rn;
    printf(" xu = %f, yu = %f\n", xu, yu);
    for (int i = 0; i < Rn; i++)
    {
        variancextmp += pow(mat.at<float>(i, 0) - xu, 2);
        varianceytmp += pow(mat.at<float>(i, 1) - yu, 2);
    }
    printf(" varixtmp = %f, variytmp = %f\n", variancextmp, varianceytmp);
    variancex = variancextmp / Rn;
    variancey = varianceytmp / Rn;
    double log1_2pivarixy = log(1 / (2 * Pi * pow(variancex, 0.5) * pow(variancey, 0.5)));
    double logL = -1 / 2 * (variancextmp / variancex + varianceytmp / variancey) + Rn * log1_2pivarixy;
    double BICparent = -2 * logL + q * log(Rn);
    printf(" logL = %lf, BICparent = %lf\n", logL, BICparent);

    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    Mat centers(Rn, 2, mat.type());
    int attempts = 1;
    TermCriteria criteria{TermCriteria::COUNT, 1, 100};
    Mat _2meanslabels(Rn, 1, labels.type());
    kmeans(mat, 2, _2meanslabels, criteria, attempts, KMEANS_PP_CENTERS, centers);

    int Rn1 = 0;
    int Rn2 = 0;
    double xu1, yu1, variancex1, variancey1, variancextmp1, varianceytmp1, xu2, yu2, variancex2, variancey2, variancextmp2, varianceytmp2 = 0.0;
    for (int i = 0; i < Rn; i++)
    {
        if (_2meanslabels.at<int>(i) == 0)
        {
            xu1 += mat.at<float>(i, 0);
            yu1 += mat.at<float>(i, 1);
            Rn1++;
        }
        else
        {
            xu2 += mat.at<float>(i, 0);
            yu2 += mat.at<float>(i, 1);
            Rn2++;
        }
    }
    xu1 = xu1 / Rn1;
    yu1 = yu1 / Rn1;
    xu2 = xu2 / Rn2;
    yu2 = yu2 / Rn2;
    printf(" Rn1 = %d, Rn2 = %d, xu1 = %lf, yu1 = %lf, xu2 = %lf, yu2 = %lf\n", Rn1, Rn2, xu1, yu1, xu2, yu2);
    //printf("xu = %f, yu = %f\n", xu, yu);
    for (int i = 0; i < Rn; i++)
    {
        if (_2meanslabels.at<int>(i) == 0)
        {
            variancextmp1 += pow(mat.at<float>(i, 0) - xu1, 2);
            varianceytmp1 += pow(mat.at<float>(i, 1) - yu1, 2);
        }
        else
        {
            variancextmp2 += pow(mat.at<float>(i, 0) - xu2, 2);
            varianceytmp2 += pow(mat.at<float>(i, 1) - yu2, 2);
        }
    }
    //printf("varixtmp = %f, variytmp = %f\n", variancextmp, varianceytmp);
    variancex1 = variancextmp1 / Rn1;
    variancey1 = varianceytmp1 / Rn1;
    variancex2 = variancextmp2 / Rn2;
    variancey2 = varianceytmp2 / Rn2;
    double log1_2pivarixy1 = log(1 / (2 * Pi * pow(variancex1, 0.5) * pow(variancey1, 0.5)));
    double log1_2pivarixy2 = log(1 / (2 * Pi * pow(variancex2, 0.5) * pow(variancey2, 0.5)));
    double logL1 = -1 / 2 * (variancextmp1 / variancex1 + varianceytmp1 / variancey1) + Rn1 * log1_2pivarixy1;
    double logL2 = -1 / 2 * (variancextmp2 / variancex2 + varianceytmp2 / variancey2) + Rn2 * log1_2pivarixy2;
    double BICchild = -2 * (logL1 + logL2) + 2 * q * log(Rn);
    printf(" logL1 = %lf, logL2 = %lf, BICchild = %lf\n", logL1, logL2, BICchild);

    if (BICparent > BICchild || (*K) > 20)
    {
        printf("split!\n");
        (*K)++;
        Mat mat1(Rn1, 2, mat.type());
        Mat mat2(Rn2, 2, mat.type());
        Mat orders1(Rn1, 1, orders.type());
        Mat orders2(Rn2, 1, orders.type());
        int Rn1tmp = 0;
        int Rn2tmp = 0;
        
        for (int i = 0; i < Rn; i++)
        {
            if (_2meanslabels.at<int>(i) == 0)
            {
                mat1.at<float>(Rn1tmp, 0) = mat.at<float>(i, 0);
                mat1.at<float>(Rn1tmp, 1) = mat.at<float>(i, 1);
                orders1.at<int>(Rn1tmp) = orders.at<int>(i);
                labels.at<int>(orders.at<int>(i)) = (*K);
                Rn1tmp++;
            }
            else
            {
                mat2.at<float>(Rn2tmp, 0) = mat.at<float>(i, 0);
                mat2.at<float>(Rn2tmp, 1) = mat.at<float>(i, 1);
                orders2.at<int>(Rn2tmp) = orders.at<int>(i);
                labels.at<int>(orders.at<int>(i)) = (*K) + 1;
                Rn2tmp++;
            }
        }
        printf("next xmeans\n");
        xmeans(mat1, orders1, labels, K);
        xmeans(mat2, orders2, labels, K);
    }
    else
    {
        return 0;
    }
}

int main(int argc, const char **argv)
{
    Mat wallimg = imread("whitebase.jpg", IMREAD_COLOR);
    Mat image_blurred_with_5x5_kernel;
    GaussianBlur(wallimg, image_blurred_with_5x5_kernel, Size(5, 5), 0);
    Mat handimg = imread("dots3.jpg", IMREAD_COLOR);
    GaussianBlur(handimg, image_blurred_with_5x5_kernel, Size(5, 5), 0);
    namedWindow("hand", WINDOW_AUTOSIZE);
    imshow("hand", handimg);

    Mat fgMaskMOG;                                  //fg mask generated by MOG method
    pMOG = bgsegm::createBackgroundSubtractorMOG(); //MOG approach
    pMOG->apply(wallimg, fgMaskMOG);
    pMOG->apply(handimg, fgMaskMOG);
    namedWindow("masked_hand", WINDOW_AUTOSIZE);
    imshow("masked_hand", fgMaskMOG);

    int whitepixel = 0;
    for (int y = 0; y < fgMaskMOG.rows; y++)
    {
        for (int x = 0; x < fgMaskMOG.cols; x++)
        {
            if (fgMaskMOG.at<unsigned char>(y, x) > 200)
            {
                whitepixel++;
            }
        }
    }
    Mat samples(whitepixel, 3, CV_32F);
    Mat orders(whitepixel, 1, CV_32F);
    whitepixel = 0;
    for (int y = 0; y < fgMaskMOG.rows; y++)
    {
        for (int x = 0; x < fgMaskMOG.cols; x++)
        {

            if (fgMaskMOG.at<unsigned char>(y, x) > 200)
            {
                samples.at<float>(whitepixel, 0) = x;
                samples.at<float>(whitepixel, 1) = y;
                orders.at<int>(whitepixel, 0) = whitepixel;
                whitepixel++;
            }
        }
    }
    printf("whitepixel = %d\n", whitepixel);
    Mat labels(whitepixel, 1, CV_32F);
    for (int i = 0; i < whitepixel; i++)
    {
        labels.at<int>(i) = 0;
    }

    xmeans(samples, orders, labels, &K);
    printf("result K = %d\n", K);
    //std::cout << "size: " << centers.size() << endl;
    //std::cout << "centers: " << centers << endl;
    //std::cout << "labels: " << labels.size() << endl;

    Mat new_image(wallimg.size(), wallimg.type());
    int cluster_idx;
    int number = 0;
    int color = 1;
    for (int y = 0; y < wallimg.rows; y++)
    {
        for (int x = 0; x < wallimg.cols; x++)
        {
            if (fgMaskMOG.at<unsigned char>(y, x) > 200)
            {
                cluster_idx = labels.at<int>(number, 0);
                new_image.at<Vec3b>(y, x)[0] = ((cluster_idx + 2) * 50) % 250; //centers.at<float>(cluster_idx, 0) * color;
                new_image.at<Vec3b>(y, x)[1] = ((cluster_idx + 2) * 130) % 250; //centers.at<float>(cluster_idx, 1) * color;
                new_image.at<Vec3b>(y, x)[2] = ((cluster_idx + 2) * 30) % 250; //centers.at<float>(cluster_idx, 2) * color;
                number++;
            }
            else
            {
                new_image.at<Vec3b>(y, x)[0] = fgMaskMOG.at<unsigned char>(y, x);
                new_image.at<Vec3b>(y, x)[1] = fgMaskMOG.at<unsigned char>(y, x);
                new_image.at<Vec3b>(y, x)[2] = fgMaskMOG.at<unsigned char>(y, x);
            }
        }
    }

    namedWindow("clustered_image", WINDOW_AUTOSIZE);
    imshow("clustered_image", new_image);
    waitKey(0);
    destroyAllWindows();

    return 0;
}