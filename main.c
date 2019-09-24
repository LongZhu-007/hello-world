#include <stdio.h>
#include <stdlib.h>
void print_1DArray(double *A, int N);
void print_2DArray(double (*A)[], int N);
void gauss_eliminate(double (*A)[], double *b, int N);
double *slover(double (*A)[], double *b, int N);
double *slover_down(double (*A)[3], double *b, int N);
void trig_decomp(double (*arr)[3], double (*L)[3], double (*U)[3], int N);
void Gauss_Jodan_eliminate(double (*Arr)[], double *b, int N);

/*��������������һ��������ʵ������Ԫ��˹��ȥ��*/
int main()
{
    printf("Hello world!\n");
    // �����������������һ������ĳߴ��С
    int N = 3;
    double A[3][3] = {1, 2, 3, 2, 5, 2, 3, 1, 5};
    double b[3] = {14, 18, 20};
    printf("����ǰ����ľ���Ϊ:\n");
    print_2DArray(A, N);
    printf("����ǰ�����b����Ϊ:\n");
    print_1DArray(b, N);
    /* ======================================== */
    // ���������������о�������Ƿֽ�
//    double L[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
//    double U[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
//    trig_decomp(A, L, U, N);
//    /* Ȼ����� LUx=b ==> Ux = a; La = b*/
//    double *a = NULL;
//    a = slover_down(L, b, N);
//    printf("���������a:\n");
//    print_1DArray(a, N);
//    double *x = NULL;
//    x = slover(U, a, N);
//    printf("���������x���Ϊ:\n");
//    print_1DArray(x, N);
//    free(a);
//    free(x);
    /* ======================================== */
    Gauss_Jodan_eliminate(A, b, N);
    printf("������˹-Լ����ȥ������֮��Ľ��:\n");
    printf("��ǰ��A����Ϊ:\n");
    print_2DArray(A, N);
    printf("��ǰ��b����Ϊ:\n");
    print_1DArray(b, N);
    // ���������ǽ����������͵Ļ���
//    printf("�����������������εĻ���:\n");
//    gauss_eliminate(A, b, N);
//    printf("�������λ���֮���ϵ������Ϊ:\n");
//    print_2DArray(A, N);
//    printf("��Ӧ��b����Ϊ:\n");
//    print_1DArray(b, N);
//    // ���������Ǳ���Ҫ�����
//    double *x = NULL;
//    x = slover(A, b, N);
//    // ��������x����Ϊ:
//    printf("������ɵ�x����Ϊ:\n");
//    print_1DArray(x, N);
//    // ������������Ҫ���ͷŶ�̬���ٵ��ڴ�ռ�
//    free(x);
    return 0;
}
/*����������һ�����㷽����ĺ���-> ���������*/
double *slover(double (*A)[3], double *b, int N){
    // ����Ҫ�����������ϵ���������Ѿ���Ϊ��������֮��ľ���
    // �������ﲢû����Ӳ��������ͻ������ݵ��ж�,
    // ������������̬�����ڴ�
    double *x = (double *)malloc(N);
    double sum_val = 0.0;
    for(int d=N-1; d>=0; d--){
        // ����ȷ�����ǵ�ǰ��������һ��Ԫ��
        sum_val = 0;
        for(int m=N-1; m> d; m--){
            sum_val += A[d][m] * x[m];
        }
        x[d] = (b[d] - sum_val) / A[d][d];
    }
    return x;
}

/* ����������һ�����㷽����ĺ���-> ���������*/
double *slover_down(double (*A)[3], double *b, int N){
    // ��̬����һ���ڴ�ռ�
    double *x = (double *)malloc(3);
    double sum_val = 0.0;
    // �������ȸ�ֵ��һ��Ԫ��
    for(int m=0; m<N; m++){ // m��ʾ��Ӧ��A�е���,
        sum_val = 0.0;
        for(int n=0; n<m; n++){ //n��ʾ��Ӧ��A�е���
            sum_val += A[m][n]*x[n];
        }
        // Ȼ���������������Ӧ��x
        x[m] = (b[m] - sum_val) / A[m][m];
    }
    return x;
}

void gauss_eliminate(double (*A)[3], double *b, int N){
    /*����Ԫ��ȥ������Ҫ������Ӧ�������ҵ���Ӧ�����ֵ��Ȼ���任����Ԫ��ȥ*/
    int max_ele_index = 0;;
    double eff = 0.0;
    for(int i=0; i<N; i++){
        max_ele_index = i;
        for(int j=i+1; j<N; j++){
            if(abs(A[j][i]) > abs(A[max_ele_index][i])){
                max_ele_index = j;
            }
        }
        /*������ѭ���󣬱�˵�������Ѿ��ҵ��˵�ǰ���ϵ����Ԫ��ֵ*/
        /*���佻������Ԫ��ȥ*/
        /*��ά����Ĵ洢��ʽ�ǰ������е�*/
        double temp = 0.0;
        for(int d=0; d<N; d++){
            temp = A[i][d];
            A[i][d] = A[max_ele_index][d];
            A[max_ele_index][d] = temp;
        }
        // ����Ҫע������ͬ����Ҫ����b�е�������Ӧλ���ϵ�ֵ
        double tet = 0.0;
        tet = *(b+i);
        *(b+i) = *(b+max_ele_index);
        *(b+max_ele_index) = tet;
        // �������֮��˵����ǰ����Ԫ�Ѿ�ѡ����
        // ������������Ҫ�����������͵Ļ���
        printf("��%d�ν���:\n", i);
        print_2DArray(A, N);
        printf("\n");
        print_1DArray(b, N);
        for(int d=i+1; d<N; d++){
            // ����Ҫ��i���ϵ�Ԫ�����õ�d����ȥ
            // �еĴ���Ӧ���Ǵ�i�п�ʼ�ĵط�
            eff = A[d][i] / A[i][i];
            for(int m=i; m<N; m++){
                A[d][m] -= A[i][m] * eff;
            }
            // ��������ͬ����Ҫ����b��Ӧλ���ϼ���
            b[d] -= eff * b[i];
        }
    }
}

/*����һ������������ӡ�����λ����*/
void print_2DArray(double (*A)[3], int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%lf\t", A[i][j]);
        }
        printf("\n");
    }
}

/*����һ������������ӡ���һά����*/
void print_1DArray(double *A, int N){
    for(int i=0; i<N; i++){
        printf("%lf\t", A[i]);
    }
    printf("\n");
}

// ������������������ϵ�к������ھ�������Ƿֽ�
// ���Ƿֽ�; triangle decomposition
void trig_decomp(double (*arr)[3], double (*L)[3], double (*U)[3], int N){
    // �����������������Ƿֽⷨ��������Ϊ��λ��������Ϊ�ǵ�λ��
    for(int i=0; i<N; i++){ // ��L�ĶԽ���Ԫ�ظ�ֵΪ1�� ��U�ĵ�һ�и�ֵΪarr�еĵ�һ��
        L[i][i] = 1;
        U[0][i] = arr[0][i];
    }
    // Ȼ�����������ȷ��L�ĵ�һ��
    for(int i=1; i<N; i++){
        // ��һ���е�����Ԫ����A�еĵ�һ�г��ϵ�һ�еĵ�һ��Ԫ��
        L[i][0] = arr[i][0] / arr[0][0];
    }
    printf("������ϳ�ʼ�ĸ�ֵ֮��Ľ��:\n");
    printf("ԭϵ������:\n");
    print_2DArray(arr, N);
    printf("L:\n");
    print_2DArray(L, N);
    printf("U:\n");
    print_2DArray(U, N);
    double sum_val = 0.0;
    // Ȼ���������ת��U��һ��, ��ʼ�������ǵ��ڵ�ǰ�е�
    // ������ѭ������U������
    for(int i=1; i<N; i++){ // �����i �Ƕ�Ӧ��U������
        printf("i=%d\n", i);
        for(int j=i; j<N; j++){ // j��ʾ����U��Ӧ������
            // ��Ӧ�����Ϊi�������������
            sum_val = 0.0;
            for(int d=i-1; d>=0; d--){ // d����Ҫ��͵�U������������������У������L[i][i] == 1,�Ӷ������U[i][j]
                sum_val += U[d][j]*L[i][d];
            }
            // ����ѭ���󣬶�Ӧ��U[i][j]���ϵ�Ԫ��,
            U[i][j] = arr[i][j] - sum_val;
        }
        printf("����U����ֵ֮��:\n");
        print_2DArray(U, N);
        // ����ѭ����˵��U�еĵ�ǰ��������
        // �����������Ǹ�ȥ���L�еĶ�Ӧ������
        // ��ʱ��i��תΪL������
        for(int d=i+1; d<N; d++){ // �˴���d��ʾ����L���У�����d==i���ϵ�Ԫ�أ���ʵ�Ѿ���ֵΪ1��
            sum_val = 0.0; // When i= 1, d=2.3.4, m=0;
            for(int m=i-1; m>=0; m--){ // �����m��������������Ҫ����ǰָ����iǰ������е����ϵ�Ԫ�����
                // L��Ӧ������U��Ӧ�������
                sum_val += L[d][m]*U[m][i];
            }
            L[d][i] = (arr[d][i] - sum_val) / U[i][i];
        }
        printf("����L����ֵ֮��:\n");
        print_2DArray(L, N);
    }
}

/* ��������ô������һ�����������ø�˹-Լ����ȥ������ʵ�ֶ��ڷ���������*/
void Gauss_Jodan_eliminate(double (*Arr)[3], double *b, int N){
    // ����ĸ�˹-Լ����ȥ�����Ƕ�������Ԫ��˹��ȥ������ǿ��ʽ
    int max_ele_index = 0;
    double effed = 0.0;
    for(int row=0; row<N; row++){ //�����ѭ������row�Ƕ�Ӧ��Arr����
        max_ele_index =row;
        for(int row_j=row+1; row_j<N; row_j++){
            if(abs(Arr[row_j][row]) >abs(Arr[max_ele_index][row])){
                max_ele_index = row_j;
            }
        }
        // ������ѭ���󣬱�õ��˶�Ӧ�������ֵ���ڵ��е�����
        // Ȼ��Ҫ����Ӧ�����ϵ�����Ԫ�ؽ�����������ȥ
        double temp =0.0;
        for(int d=row; d<N; d++){
            temp = Arr[max_ele_index][d];
            Arr[max_ele_index][d] = Arr[row][d];
            Arr[row][d] = temp;
        }
        // ����ԪԪ�ؽ�����Ϻ�������Ҫ��b�ж�Ӧλ���ϵ�Ԫ�ؽ��н���
        temp = b[max_ele_index];
        b[max_ele_index] = b[row];
        b[row] = temp;
        printf("��%d���޸Ĺ����Arr����:\n", row);
        print_2DArray(Arr, 3);
        printf("��%d���޸Ĺ����b����:\n", row);
        print_1DArray(b, 3);
        // ������������Ҫ�����ǣ�����ÿ���ϵĻ������
        for(int r=0; r<N; r++){
            // ��ǰ��������row, ��ǰ������������r
            if(r==row){
                //������������ǵ�ǰ������,��ôϵ�����ǵ�һ��Ԫ��
                effed = Arr[row][row];
                // Ȼ�����ø�ϵ����������ǰ���е������е�Ԫ��
                for(int h=row; h<N; h++){
                    Arr[row][h] /= effed;
                }
                // Ȼ���������������Ҫ������b�ж�Ӧ����λ���ϵ�Ԫ��
                b[row] /= effed;
            }else{
                // ��ǰ��ϵ��effed����Ҫ�뵱ǰ�����н���һ���Ա���
                effed = Arr[r][row] / Arr[row][row];
                // Ȼ�����ø�ϵ����������ǰ���е������е�Ԫ��
                for(int h=row; h<N; h++){ // ������Ϊ��Ԫ���ڵ��е�row��ǰ���Ԫ�ض�Ϊ0,��������ı������Դ�row��ʼ
                    Arr[r][h] -= Arr[row][h] *effed;
                }
                // ����������Ҫ����b�ж�Ӧ����λ���ϵ�Ԫ��
                b[r] -= effed*b[row];
            }
        }
        printf("�޸�������֮���Arr����:\n");
        print_2DArray(Arr, 3);
        printf("�޸�������֮���b����:\n");
        print_1DArray(b, 3);
    }
}

