#include <stdio.h>
#include <stdlib.h>
void print_1DArray(double *A, int N);
void print_2DArray(double (*A)[], int N);
void gauss_eliminate(double (*A)[], double *b, int N);
double *slover(double (*A)[], double *b, int N);
double *slover_down(double (*A)[3], double *b, int N);
void trig_decomp(double (*arr)[3], double (*L)[3], double (*U)[3], int N);
void Gauss_Jodan_eliminate(double (*Arr)[], double *b, int N);

/*本程序是来创建一个方法，实现列主元高斯消去法*/
int main()
{
    printf("Hello world!\n");
    // 我们在这合理来创建一个矩阵的尺寸大小
    int N = 3;
    double A[3][3] = {1, 2, 3, 2, 5, 2, 3, 1, 5};
    double b[3] = {14, 18, 20};
    printf("您当前输入的矩阵为:\n");
    print_2DArray(A, N);
    printf("您当前输入的b矩阵为:\n");
    print_1DArray(b, N);
    /* ======================================== */
    // 接下来我们来进行矩阵的三角分解
//    double L[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
//    double U[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
//    trig_decomp(A, L, U, N);
//    /* 然后根据 LUx=b ==> Ux = a; La = b*/
//    double *a = NULL;
//    a = slover_down(L, b, N);
//    printf("计算出来的a:\n");
//    print_1DArray(a, N);
//    double *x = NULL;
//    x = slover(U, a, N);
//    printf("计算出来的x结果为:\n");
//    print_1DArray(x, N);
//    free(a);
//    free(x);
    /* ======================================== */
    Gauss_Jodan_eliminate(A, b, N);
    printf("经过高斯-约旦消去法处理之后的结果:\n");
    printf("当前的A矩阵为:\n");
    print_2DArray(A, N);
    printf("当前的b矩阵为:\n");
    print_1DArray(b, N);
    // 接下来我们进行上三角型的划分
//    printf("接下来进行上三角形的划分:\n");
//    gauss_eliminate(A, b, N);
//    printf("上三角形划分之后的系数矩阵为:\n");
//    print_2DArray(A, N);
//    printf("对应的b矩阵为:\n");
//    print_1DArray(b, N);
//    // 接下来我们便需要来求解
//    double *x = NULL;
//    x = slover(A, b, N);
//    // 计算过后的x向量为:
//    printf("计算完成的x矩阵为:\n");
//    print_1DArray(x, N);
//    // 接下来我们需要来释放动态开辟的内存空间
//    free(x);
    return 0;
}
/*接下来创建一个计算方程组的函数-> 上三角求解*/
double *slover(double (*A)[3], double *b, int N){
    // 这里要求您所输入的系数矩阵是已经化为上三角型之后的矩阵
    // 但是这里并没有添加参数的类型或者数据的判断,
    // 这里我们来动态开辟内存
    double *x = (double *)malloc(N);
    double sum_val = 0.0;
    for(int d=N-1; d>=0; d--){
        // 首先确定的是当前结果的最后一个元素
        sum_val = 0;
        for(int m=N-1; m> d; m--){
            sum_val += A[d][m] * x[m];
        }
        x[d] = (b[d] - sum_val) / A[d][d];
    }
    return x;
}

/* 接下来创建一个计算方程组的函数-> 下三角求解*/
double *slover_down(double (*A)[3], double *b, int N){
    // 动态开辟一个内存空间
    double *x = (double *)malloc(3);
    double sum_val = 0.0;
    // 这里首先赋值第一个元素
    for(int m=0; m<N; m++){ // m表示对应的A中的行,
        sum_val = 0.0;
        for(int n=0; n<m; n++){ //n表示对应的A中的列
            sum_val += A[m][n]*x[n];
        }
        // 然后我们来计算出对应的x
        x[m] = (b[m] - sum_val) / A[m][m];
    }
    return x;
}

void gauss_eliminate(double (*A)[3], double *b, int N){
    /*列主元消去法，需要遍历对应的列来找到对应的最大值，然后将其换到主元上去*/
    int max_ele_index = 0;;
    double eff = 0.0;
    for(int i=0; i<N; i++){
        max_ele_index = i;
        for(int j=i+1; j<N; j++){
            if(abs(A[j][i]) > abs(A[max_ele_index][i])){
                max_ele_index = j;
            }
        }
        /*待跳出循环后，便说明我们已经找到了当前列上的最大元素值*/
        /*将其交换到主元上去*/
        /*二维数组的存储方式是按行排列的*/
        double temp = 0.0;
        for(int d=0; d<N; d++){
            temp = A[i][d];
            A[i][d] = A[max_ele_index][d];
            A[max_ele_index][d] = temp;
        }
        // 这里要注意我们同样需要交换b中的两个对应位置上的值
        double tet = 0.0;
        tet = *(b+i);
        *(b+i) = *(b+max_ele_index);
        *(b+max_ele_index) = tet;
        // 交换完成之后，说明当前的主元已经选定了
        // 接下来我们需要进行上三角型的划分
        printf("第%d次交换:\n", i);
        print_2DArray(A, N);
        printf("\n");
        print_1DArray(b, N);
        for(int d=i+1; d<N; d++){
            // 我们要把i行上的元素作用到d行上去
            // 列的处理应该是从i列开始的地方
            eff = A[d][i] / A[i][i];
            for(int m=i; m<N; m++){
                A[d][m] -= A[i][m] * eff;
            }
            // 我们这里同样需要进行b对应位置上计算
            b[d] -= eff * b[i];
        }
    }
}

/*创建一个方法，来打印输出二位数组*/
void print_2DArray(double (*A)[3], int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%lf\t", A[i][j]);
        }
        printf("\n");
    }
}

/*创建一个方法，来打印输出一维数组*/
void print_1DArray(double *A, int N){
    for(int i=0; i<N; i++){
        printf("%lf\t", A[i]);
    }
    printf("\n");
}

// 接下来，我们来创建系列函数用于矩阵的三角分解
// 三角分解; triangle decomposition
void trig_decomp(double (*arr)[3], double (*L)[3], double (*U)[3], int N){
    // 我们这里来进行三角分解法，下三角为单位阵，上三角为非单位阵
    for(int i=0; i<N; i++){ // 给L的对角线元素赋值为1， 给U的第一行赋值为arr中的第一行
        L[i][i] = 1;
        U[0][i] = arr[0][i];
    }
    // 然后接下来我们确定L的第一列
    for(int i=1; i<N; i++){
        // 第一列中的所有元素是A中的第一列除上第一行的第一个元素
        L[i][0] = arr[i][0] / arr[0][0];
    }
    printf("进行完毕初始的赋值之后的结果:\n");
    printf("原系数矩阵:\n");
    print_2DArray(arr, N);
    printf("L:\n");
    print_2DArray(L, N);
    printf("U:\n");
    print_2DArray(U, N);
    double sum_val = 0.0;
    // 然后接下来便转到U下一行, 起始操作列是等于当前行的
    // 接下来循环的是U的行数
    for(int i=1; i<N; i++){ // 这里的i 是对应的U的行数
        printf("i=%d\n", i);
        for(int j=i; j<N; j++){ // j表示的是U对应的列数
            // 对应的求和为i的上面的所有行
            sum_val = 0.0;
            for(int d=i-1; d>=0; d--){ // d是需要求和的U的行数的上面的所有行，这里的L[i][i] == 1,从而计算出U[i][j]
                sum_val += U[d][j]*L[i][d];
            }
            // 跳出循环后，对应的U[i][j]列上的元素,
            U[i][j] = arr[i][j] - sum_val;
        }
        printf("对于U赋完值之后:\n");
        print_2DArray(U, N);
        // 跳出循环后，说明U中的当前行填充完毕
        // 接下来，我们该去填充L中的对应的列了
        // 此时的i将转为L的列数
        for(int d=i+1; d<N; d++){ // 此处的d表示的是L的行，对于d==i行上的元素，其实已经赋值为1了
            sum_val = 0.0; // When i= 1, d=2.3.4, m=0;
            for(int m=i-1; m>=0; m--){ // 这里的m是列数，我们需要将当前指定的i前面的所有的列上的元素求和
                // L对应的行与U对应的列相乘
                sum_val += L[d][m]*U[m][i];
            }
            L[d][i] = (arr[d][i] - sum_val) / U[i][i];
        }
        printf("对于L赋完值之后:\n");
        print_2DArray(L, N);
    }
}

/* 在这里我么来创建一个函数，采用高斯-约旦消去法，来实现对于方程组的求解*/
void Gauss_Jodan_eliminate(double (*Arr)[3], double *b, int N){
    // 这里的高斯-约旦消去法，是对于列主元高斯消去法的增强形式
    int max_ele_index = 0;
    double effed = 0.0;
    for(int row=0; row<N; row++){ //这里的循环变量row是对应的Arr的行
        max_ele_index =row;
        for(int row_j=row+1; row_j<N; row_j++){
            if(abs(Arr[row_j][row]) >abs(Arr[max_ele_index][row])){
                max_ele_index = row_j;
            }
        }
        // 待跳出循环后，便得到了对应的最大主值所在的行的索引
        // 然后要将对应的行上的所有元素交换到主行上去
        double temp =0.0;
        for(int d=row; d<N; d++){
            temp = Arr[max_ele_index][d];
            Arr[max_ele_index][d] = Arr[row][d];
            Arr[row][d] = temp;
        }
        // 待主元元素交换完毕后，我们需要把b中对应位置上的元素进行交换
        temp = b[max_ele_index];
        b[max_ele_index] = b[row];
        b[row] = temp;
        printf("第%d次修改过后的Arr矩阵:\n", row);
        print_2DArray(Arr, 3);
        printf("第%d次修改过后的b矩阵:\n", row);
        print_1DArray(b, 3);
        // 接下来，我们要做的是，进行每行上的化简操作
        for(int r=0; r<N; r++){
            // 当前的主行是row, 当前遍历到的行是r
            if(r==row){
                //如果遍历到的是当前的主行,那么系数便是第一个元素
                effed = Arr[row][row];
                // 然后利用该系数，来处理当前行中的其余列的元素
                for(int h=row; h<N; h++){
                    Arr[row][h] /= effed;
                }
                // 然后接下来，我们需要来处理b中对应索引位置上的元素
                b[row] /= effed;
            }else{
                // 当前的系数effed则需要与当前的主行进行一个对比了
                effed = Arr[r][row] / Arr[row][row];
                // 然后利用该系数，来处理当前行中的其余列的元素
                for(int h=row; h<N; h++){ // 这里因为主元所在的行的row列前面的元素都为0,所以这里的遍历可以从row开始
                    Arr[r][h] -= Arr[row][h] *effed;
                }
                // 接下来，需要处理b中对应索引位置上的元素
                b[r] -= effed*b[row];
            }
        }
        printf("修改其他行之后的Arr矩阵:\n");
        print_2DArray(Arr, 3);
        printf("修改其他行之后的b矩阵:\n");
        print_1DArray(b, 3);
    }
}

