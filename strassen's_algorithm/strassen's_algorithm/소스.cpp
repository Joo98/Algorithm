#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#define MAX_SIZE 10
using namespace std;
//행렬의 최대 크기는 10으로 한다

class matrix
{
private:
    int a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE][MAX_SIZE], result[MAX_SIZE][MAX_SIZE];
    int n;
    long int start_time, end_time;
    void strassen(int n, int a1[][MAX_SIZE], int b1[][MAX_SIZE], int                                             c1[][MAX_SIZE]);
public:
    void init();

    const void print(const char* msg);


    void multiplication();
    void strassen_start();
};



//쉬트라센 알고리즘 시작
void matrix::strassen_start()
{
    double start, stop, timer;
    start = clock();
    //수행 시간을 구하기 위해 start에 시작 시간을 넣는다  
    int i, j;



    //계산된 결과가 저장될 배열을 만든다
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = 0;
    strassen(n, a, b, result);
    //쉬트라센 함수를 호출한다
    stop = clock();
    timer = ((double)(stop - start)) / CLK_TCK;
    cout << timer << endl;
    //쉬트라센 함수를 호출후 계산을 마치고 돌아온 시간을 기록하여 수행시간을 계산한다
}





//쉬트라센 함수
void matrix::strassen(int n, int a1[][MAX_SIZE], int b1[][MAX_SIZE], int c1[][MAX_SIZE])
{

    int i, j;
    int m1[MAX_SIZE][MAX_SIZE], m2[MAX_SIZE][MAX_SIZE],

        m3[MAX_SIZE][MAX_SIZE], m4[MAX_SIZE][MAX_SIZE],

        m5[MAX_SIZE][MAX_SIZE], m6[MAX_SIZE][MAX_SIZE],

        m7[MAX_SIZE][MAX_SIZE];


    int temp1[MAX_SIZE][MAX_SIZE], temp2[MAX_SIZE][MAX_SIZE];

    int sm1, sm2, sm3, sm4, sm5, sm6, sm7;

    if (n <= 2) {
        //만약 행렬의 크기가 2라면 바로 곱셈을 한다
        //그 이유는 손익분기점이 2이기 때문이다 
        sm1 = (a1[0][0] + a1[1][1]) * (b1[0][0] + b1[1][1]);
        sm2 = (a1[1][0] + a1[1][1]) * b1[0][0];
        sm3 = a1[0][0] * (b1[0][1] - b1[1][1]);
        sm4 = a1[1][1] * (b1[1][0] - b1[0][0]);
        sm5 = (a1[0][0] + a1[0][1]) * b1[1][1];
        sm6 = (a1[1][0] - a1[0][0]) * (b1[0][0] + b1[0][1]);
        sm7 = (a1[0][1] - a1[1][1]) * (b1[1][0] + b1[1][1]);

        c1[0][0] = sm1 + sm4 - sm5 + sm7;
        c1[0][1] = sm3 + sm5;
        c1[1][0] = sm2 + sm4;
        c1[1][1] = sm1 + sm3 - sm2 + sm6;
    }
    else {
        //손익 분기점을 넘는다면    
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i][j] + a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j] + b1[i + n / 2][j + n / 2];
        strassen(n / 2, temp1, temp2, m1);
        //을 계산하여 쉬트라센 호출
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i + n / 2][j] + a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j];
        strassen(n / 2, temp1, temp2, m2);
        //계산하여 쉬트라센호출
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i][j];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j + n / 2] - b1[i + n / 2][j + n / 2];
        strassen(n / 2, temp1, temp2, m3);
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i + n / 2][j] - b1[i][j];
        strassen(n / 2, temp1, temp2, m4);
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i][j] + a1[i][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i + n / 2][j + n / 2];
        strassen(n / 2, temp1, temp2, m5);
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i + n / 2][j] - a1[i][j];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j] + b1[i][j + n / 2];
        strassen(n / 2, temp1, temp2, m6);
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i][j + n / 2] - a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i + n / 2][j] + b1[i + n / 2][j + n / 2];
        strassen(n / 2, temp1, temp2, m7);
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                c1[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                c1[i][j + n / 2] = m3[i][j] + m5[i][j];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                c1[i + n / 2][j] = m2[i][j] + m4[i][j];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                c1[i + n / 2][j + n / 2] = m1[i][j] + m3[i][j] - m2[i][j] + m6[i][j];
        //이런 식으로 계산하여 각각 쉬트라센을 호출한다
    }
}

void matrix::init()
//행렬의 크기를 결정
{
    int i, j;
    do {
        cout << endl << endl;
        cout << "행렬의 크기를 결정 (2의 제곱수)" << endl;
        cout << ">> ";
        cin >> n;
        //행렬의 크기를 결정하기 위해 사용자로부터 행열의 크기를 부여받는다
        if (n > MAX_SIZE) continue;
        i = n;
        j = 0;
        do {
            i = i / 2;
        } while (i != 1);

        if (j == 0) break;
    } while (1);
    srand((unsigned)time(NULL));
    //시간함수를 호출하여 랜덤함수를 호출한다


    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            a[i][j] = rand() % 10;
            b[i][j] = rand() % 10;
            //각각의 랜덤함수는 10보다 작은수로 결정한다.
        }

    cout << "================================" << endl;
    printf("A\n\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d  ", a[i][j]);
        }
        printf("\n\n");
    }
    cout << "================================" << endl;
    printf("B\n\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d  ", b[i][j]);
        }
        printf("\n\n");
    }
    
}
void matrix::multiplication()
//표준 알고리즘을 이용하여 행렬을 계산한다
{
    double start, stop, timer;
    start = clock();
    //역시 시간을 계산하기위해 start에 현재 시간을 넣고
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            result[i][j] = 0;
            for (k = 0; k < n; k++)
                result[i][j] = result[i][j] + a[i][k] * b[k][j];
        }
    //각각 행렬을 가로세로로 각각 곱해서 표준 알고리즘으로 계산
    stop = clock();
    timer = ((double)(stop - start)) / CLK_TCK;
    cout << timer << endl;
    //걸린 시간을 출력한다
}

const void matrix::print(const char* msg)
{
    int i, j;

    cout << "============" << endl << msg << endl << "============" << endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cout.width(7);
            cout << result[i][j];
        }
        cout << endl;
    }


}


int main()
//main 이지만 하는일은 모든 함수들의 계산한 결과값을 출력하는 것만 수행
{
    matrix m;

    m.init();
    m.multiplication();
    m.print("표준 알고리즘");
    m.strassen_start();
    m.print("쉬트라센 알고리즘");
    return 0;
}