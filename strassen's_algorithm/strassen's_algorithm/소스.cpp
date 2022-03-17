#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#define MAX_SIZE 10
using namespace std;
//����� �ִ� ũ��� 10���� �Ѵ�

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



//��Ʈ�� �˰��� ����
void matrix::strassen_start()
{
    double start, stop, timer;
    start = clock();
    //���� �ð��� ���ϱ� ���� start�� ���� �ð��� �ִ´�  
    int i, j;



    //���� ����� ����� �迭�� �����
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = 0;
    strassen(n, a, b, result);
    //��Ʈ�� �Լ��� ȣ���Ѵ�
    stop = clock();
    timer = ((double)(stop - start)) / CLK_TCK;
    cout << timer << endl;
    //��Ʈ�� �Լ��� ȣ���� ����� ��ġ�� ���ƿ� �ð��� ����Ͽ� ����ð��� ����Ѵ�
}





//��Ʈ�� �Լ�
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
        //���� ����� ũ�Ⱑ 2��� �ٷ� ������ �Ѵ�
        //�� ������ ���ͺб����� 2�̱� �����̴� 
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
        //���� �б����� �Ѵ´ٸ�    
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i][j] + a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j] + b1[i + n / 2][j + n / 2];
        strassen(n / 2, temp1, temp2, m1);
        //�� ����Ͽ� ��Ʈ�� ȣ��
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp1[i][j] = a1[i + n / 2][j] + a1[i + n / 2][j + n / 2];
        for (i = 0; i < n / 2; i++)
            for (j = 0; j < n / 2; j++)
                temp2[i][j] = b1[i][j];
        strassen(n / 2, temp1, temp2, m2);
        //����Ͽ� ��Ʈ��ȣ��
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
        //�̷� ������ ����Ͽ� ���� ��Ʈ���� ȣ���Ѵ�
    }
}

void matrix::init()
//����� ũ�⸦ ����
{
    int i, j;
    do {
        cout << endl << endl;
        cout << "����� ũ�⸦ ���� (2�� ������)" << endl;
        cout << ">> ";
        cin >> n;
        //����� ũ�⸦ �����ϱ� ���� ����ڷκ��� �࿭�� ũ�⸦ �ο��޴´�
        if (n > MAX_SIZE) continue;
        i = n;
        j = 0;
        do {
            i = i / 2;
        } while (i != 1);

        if (j == 0) break;
    } while (1);
    srand((unsigned)time(NULL));
    //�ð��Լ��� ȣ���Ͽ� �����Լ��� ȣ���Ѵ�


    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            a[i][j] = rand() % 10;
            b[i][j] = rand() % 10;
            //������ �����Լ��� 10���� �������� �����Ѵ�.
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
//ǥ�� �˰����� �̿��Ͽ� ����� ����Ѵ�
{
    double start, stop, timer;
    start = clock();
    //���� �ð��� ����ϱ����� start�� ���� �ð��� �ְ�
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
    //���� ����� ���μ��η� ���� ���ؼ� ǥ�� �˰������� ���
    stop = clock();
    timer = ((double)(stop - start)) / CLK_TCK;
    cout << timer << endl;
    //�ɸ� �ð��� ����Ѵ�
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
//main ������ �ϴ����� ��� �Լ����� ����� ������� ����ϴ� �͸� ����
{
    matrix m;

    m.init();
    m.multiplication();
    m.print("ǥ�� �˰���");
    m.strassen_start();
    m.print("��Ʈ�� �˰���");
    return 0;
}