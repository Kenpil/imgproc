#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/bgsegm.hpp>
#include <opencv2/opencv.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/format.hpp>
using namespace cv;
using namespace std;

Ptr<BackgroundSubtractor> pMOG; //MOG Background subtractor
//コンパイル: g++ test.cpp -I/usr/local/include/opencv2 -I/usr/local/include/opencv -L/usr/lib -lopencv_core -lopencv_imgcodecs -lopencv_highgui -lopencv_imgproc -lopencv_bgsegm

void show_result(const Mat &labels, const Mat &centers, int cluster_numger, int height, int width);

int main(int argc, const char **argv)
{
    Mat wallimg = imread("wall.jpeg", IMREAD_COLOR);
    Mat image_blurred_with_5x5_kernel;
    GaussianBlur(wallimg, image_blurred_with_5x5_kernel, Size(5, 5), 0);
    Mat handimg = imread("hand.jpeg", IMREAD_COLOR);
    GaussianBlur(handimg, image_blurred_with_5x5_kernel, Size(5, 5), 0);
    //namedWindow("wall", WINDOW_AUTOSIZE);
    //imshow("wall", wallimg);
    namedWindow("hand", WINDOW_AUTOSIZE);
    imshow("hand", handimg);

    Mat fgMaskMOG;                                  //fg mask generated by MOG method
    pMOG = bgsegm::createBackgroundSubtractorMOG(); //MOG approach
    pMOG->apply(wallimg, fgMaskMOG);
    pMOG->apply(handimg, fgMaskMOG);
    namedWindow("masked_hand", WINDOW_AUTOSIZE);
    imshow("masked_hand", fgMaskMOG);

    //fgMaskMOG = wallimg;
    //printf("rows = %d, cols = %d\n", fgMaskMOG.rows, fgMaskMOG.cols);
    //printf("%d %d %d\n", fgMaskMOG.at<Vec3b>(0, 280)[0], fgMaskMOG.at<Vec3b>(0, 280)[1], fgMaskMOG.at<Vec3b>(0, 280)[2]);
    Mat samples(fgMaskMOG.rows * fgMaskMOG.cols, 3, CV_32F);
    //Mat samples = fgMaskMOG;
    int calcsize = 0;
    for (int y = 0; y < fgMaskMOG.rows; y++)
    {
        for (int x = 0; x < fgMaskMOG.cols; x++)
        {
            //for (int y = 0; y < wallimg.rows; y++)
            // {
            //for (int x = 0; x < wallimg.cols; x++)
            //{
            for (int z = 0; z < 3; z++)
            {
                if (fgMaskMOG.at<unsigned char>(y, x) > 200)
                {
                    samples.at<float>(y + x * fgMaskMOG.rows, z) = fgMaskMOG.at<unsigned char>(y, x);
                    //samples.at<float>(calcsize, z) = fgMaskMOG.at<unsigned char>(y, x);
                    calcsize++;
                }else{
                    samples.at<float>(y + x * fgMaskMOG.rows, z) = 0;
                }
            }
        }
    }
    printf("calcsize = %d\n", calcsize);

    int clusterCount = 3;
    Mat labels;
    int attempts = 1;
    Mat centers(clusterCount, 1, fgMaskMOG.type());
    TermCriteria criteria{TermCriteria::COUNT, 1, 100};

    //kmeans(samples, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers);
    kmeans(samples, clusterCount, labels, criteria, attempts, KMEANS_PP_CENTERS, centers);

    // Mat new_image(fgMaskMOG.size(), fgMaskMOG.type());
    Mat new_image(wallimg.size(), wallimg.type());
    int cluster_idx;
    //printf("rows = %d, cols = %d\n", fgMaskMOG.rows, fgMaskMOG.cols);
    //printf("%f %f %f\n", centers.at<float>(cluster_idx, 0), centers.at<float>(cluster_idx, 1), centers.at<float>(cluster_idx, 2));

    //for (int y = 0; y < fgMaskMOG.rows; y++)
    //{
    //for (int x = 0; x < fgMaskMOG.cols; x++)
    //{
    for (int y = 0; y < wallimg.rows; y++)
    {
        for (int x = 0; x < wallimg.cols; x++)
        {

            cluster_idx = labels.at<int>(y + x * fgMaskMOG.rows, 0);
            new_image.at<Vec3b>(y, x)[0] = centers.at<float>(cluster_idx * 0.5, 0);
            new_image.at<Vec3b>(y, x)[1] = centers.at<float>(cluster_idx * 0.1, 1);
            new_image.at<Vec3b>(y, x)[2] = centers.at<float>(cluster_idx, 2);
            if (y == 0)
            {
                printf("idx = %d, %d  ", cluster_idx, x);
            }
        }
        //printf("idx = %d", cluster_idx);
    }

    namedWindow("clustered_image", WINDOW_AUTOSIZE);
    imshow("clustered_image", new_image);
    //*/
    waitKey(0);
    destroyAllWindows();

    return 0;
}

void show_result(const Mat &labels, const Mat &centers, int cluster_number, int height, int width)
{
    /*
     画像の形状を(H,W)とすると、labelsの形状は(H*W,1)である。各行にはラベルの値が整数値で収められている。
     cluster_numberを3とすれば、0,1,2のいずれかとなる。
     centersの形状は(cluster_numeber,3)である。各行にクラスタの中心座標が収められている。今の場合、(B,G,R)の
     値が並ぶことになる。
     */
    std::cout << "labels: "
              << "(" << labels.rows << ", " << labels.cols << ")" << std::endl;
    std::cout << "centers: "
              << "(" << centers.rows << ", " << centers.cols << ")" << std::endl;
    assert(labels.type() == CV_32SC1);
    assert(centers.type() == CV_32FC1);

    /*
     centersの要素の型は浮動小数点である。これを0から255までのstd::uint8_tの型に変更する。
     さらに、チャンネル数を3に変更する。
     */
    Mat centers_u8{};
    centers.convertTo(centers_u8, CV_8UC1, 255.0);
    const auto centers_u8c3 = centers_u8.reshape(3);
    //assert(centers_u8c3.type() == CV_8UC3);

    // K-meansクラスタリングの結果を画像に変換して表示するので、入れ物を初期化する。
    //Mat rgb_image(height, width, CV_8UC3);

    /*
     画像に変換する。
        1. labelsからひとつずつ値を取り出す。
        2. そのラベルに相当するRGB値を取り出して、rgb_imageへコピーする。
     cv::Matの要素へのアクセスは下記のようにポインタを使うのが最速である。
     */
    //boost::copy(
    //    boost::make_iterator_range(labels.begin<int>(), labels.end<int>()) | boost::adaptors::transformed([&centers_u8c3](const auto &label) { return *centers_u8c3.ptr<Vec3b>(label); }),
    //    rgb_image.begin<Vec3b>());

    // 表示する。
    //imshow("result", rgb_image);

    // 出力パスを作り、保存する。
    //const auto output_path = (boost::format("%1%_result.jpg") % cluster_number).str();
    //imwrite(output_path, rgb_image);

    //waitKey();
}
