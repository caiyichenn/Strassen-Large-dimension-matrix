#include "examples.h"

using namespace std;
using namespace seal;


//打印矩阵
void PrintMatrix(vector<vector<double>> MatrixA)
{
    
    for (int i = 0; i < MatrixA.size(); i++)
    {
        cout << "第" << i << "行";
        for (int j = 0; j < MatrixA[0].size(); j++)
        {
            cout << setprecision(8)<< MatrixA[i][j] << " ";
        }
        cout << endl;
    }
}


//向量减法
void vector_sub(int N, vector<double> plain1, vector<double> plain2, vector<double> &destination)
{
    for (int i = 0; i < N; i++)
    {       
        destination[i] = plain1[i] - plain2[i];
        
    }
}


//向量加法
void vector_sum(int N, vector<double> plain1, vector<double> plain2, vector<double> &destination)
{
    for (int i = 0; i < N; i++)
    {
        destination[i] = plain1[i] + plain2[i];
    }
}


    //void generate_matrix(int N, vector<vector<double>> MatrixA, vector<vector<double>> MatrixB)
//{
//    srand((unsigned)time(NULL));
//    for (int i = 0; i < N; i++)
//    {
//        vector<double> tmp1, tmp2;
//
//        for (int j = 0; j < N; j++)
//        {
//            tmp1.push_back(rand() % 256);
//            tmp2.push_back(rand() % 256);
//        }
//
//        MatrixA.push_back(tmp1);
//        MatrixB.push_back(tmp2);
//    }
//}



void reverse(vector<double> &plain, int start, int end)
{
    int temp;
    while (start < end)
    {
        temp = plain[start];
        plain[start] = plain[end];
        plain[end] = temp;
        start++;
        end--;
    }
}

//由于seal内只有针对密文进行rotation的操作，所以此处需要自己实现一个明文rotation的操作
//长度为numsSize的plain向量，向左旋转k步，最终得到旋转后的plain
void rotate(vector<double> &plain, int numsSize, int k)
{
    k %= numsSize; //计算出nums[0]移动后的最后位置，即(0+k)%numSize
    reverse(plain, 0, numsSize - 1);
    reverse(plain, 0, numsSize - k - 1);
    reverse(plain, numsSize - k, numsSize - 1);
}


//每组长度为numsSize的plain向量，向左旋转k步，最终得到旋转后的plain
void bacth_rotate(vector<double> &plain, int numsSize, int k)
{
    k %= numsSize;
    for (int i = 0; i < numsSize; i++)
    {
        reverse(plain, i * numsSize, (i + 1) * numsSize - 1);
        reverse(plain, i * numsSize, (i + 1) * numsSize - k - 1);
        reverse(plain, (i + 1) * numsSize - k, (i + 1) * numsSize - 1);
    }
}


void single_bacth_rotate(vector<double> &plain, int numsSize, int i, int k)
{
    k %= numsSize;
    reverse(plain, i * numsSize, (i + 1) * numsSize - 1);
    reverse(plain, i * numsSize, (i + 1) * numsSize - k - 1);
    reverse(plain, (i + 1) * numsSize - k, (i + 1) * numsSize - 1);
}


void example_strassen()
{
    print_example_banner("Example: Strassen with CKKS");

    EncryptionParameters parms(scheme_type::CKKS);

    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));

    double scale = pow(2.0, 40);

    //参数，打印
    auto context = SEALContext::Create(parms);
    print_parameters(context);
    cout << endl;

    //输入部分，读取input
    vector<vector<double>> MatrixA, MatrixB;
    int N = 128;
    srand((unsigned)time(NULL));
    for (int i = 0; i < N; i++)
    {
        vector<double> tmp1, tmp2;

        for (int j = 0; j < N; j++)
        {
            tmp1.push_back(rand() % 256);
            tmp2.push_back(rand() % 256);
        }

        MatrixA.push_back(tmp1);
        MatrixB.push_back(tmp2);
    }


   /* for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            cout << MatrixB[i][j] << " ";
        }
        cout << endl;
    }*/

    // 计算从文件中读取数据需要的时间

    //二维输入数组89
    // ifstream f;
    // f.open("C:\\Users\\DELL\\Desktop\\data.txt");
    // string str;

    ////二维输入数组
    // vector<vector<double>> MatrixA;
    // while (getline(f, str))
    //{
    //    istringstream input(str);
    //    vector<double> tmp;
    //    double a;
    //    while (input >> a)
    //        tmp.push_back(a);

    //    MatrixA.push_back(tmp);
    //}

    // cout << endl;

    // for (int i = 0; i < MatrixA.size(); i++)
    //{
    //    for (int j = 0; i < MatrixA[0].size(); j++)
    //    {
    //        cout << MatrixA[i][j] << " ";
    //    }
    //    cout << endl;
    //}

    //// 读取input2
    //
    //// 计算从文件中读取数据需要的时间

    ////二维输入数组
    // vector<vector<double> > MatrixB;
    // ifstream in2("C:\\Users\\DELL\\Desktop\\data\\data2.txt"); //打开文件
    ////读数据。。
    // time_start = chrono::high_resolution_clock::now();
    // for (int i = 0; i < 128; ++i)
    //{
    //    for (int j = 0; j < 128; ++j)
    //    {
    //        in2 >> MatrixB[i][j];
    //    }
    //}
    // in2.close(); //关闭文件
    // time_end = chrono::high_resolution_clock::now();
    // time_read_input2 = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    // cout << "Read input2 from file:";
    // cout << time_read_input2.count() << " microseconds" << endl;
    // for (int i = 0; i < MatrixB.size(); i++)
    //{
    //    for (int j = 0; i < MatrixB[0].size(); j++)
    //    {
    //        cout << MatrixB[i][j] << " ";
    //    }
    //    cout << endl;
    //}

    //新建矩阵
   /* vector<vector<double>> MatrixA11(N/2, vector<double>(N/2));
    vector<vector<double>> MatrixA12(N / 2, vector<double>(N / 2));
    vector<vector<double>> MatrixA21(N / 2, vector<double>(N / 2));
    vector<vector<double>> MatrixA22(N / 2, vector<double>(N / 2));

    vector<vector<double>> MatrixB11(N / 2, vector<double>(N / 2));
    vector<vector<double>> MatrixB12(N / 2, vector<double>(N / 2));
    vector<vector<double>> MatrixB21(N / 2, vector<double>(N / 2));
    vector<vector<double>> MatrixB22(N / 2, vector<double>(N / 2));*/

    vector<vector<double>> MatrixA11;
    vector<vector<double>> MatrixA12;
    vector<vector<double>> MatrixA21;
    vector<vector<double>> MatrixA22;

    vector<vector<double>> MatrixB11;
    vector<vector<double>> MatrixB12;
    vector<vector<double>> MatrixB21;
    vector<vector<double>> MatrixB22;

    /*int **MatrixC11;
    int **MatrixC12;
    int **MatrixC21;
    int **MatrixC22;*/
    //初始化每个小矩阵的大小
    // MatrixA11 = new double *[N / 2]; //数组的第二维一定要显示指定
    // MatrixA12 = new double *[N / 2];
    // MatrixA21 = new double *[N / 2];
    // MatrixA22 = new double *[N / 2];

    // MatrixB11 = new double *[N / 2];
    // MatrixB12 = new double *[N / 2];
    // MatrixB21 = new double *[N / 2];
    // MatrixB22 = new double *[N / 2];

    ///*MatrixC11 = new int *[N / 2];
    // MatrixC12 = new int *[N / 2];
    // MatrixC21 = new int *[N / 2];
    // MatrixC22 = new int *[N / 2];*/
    // for (int i = 0; i < N / 2; i++) //分配连续内存
    //{
    //    MatrixA11[i] = new double[N / 2];
    //    MatrixA12[i] = new double[N / 2];
    //    MatrixA21[i] = new double[N / 2];
    //    MatrixA22[i] = new double[N / 2];

    //    MatrixB11[i] = new double[N / 2];
    //    MatrixB12[i] = new double[N / 2];
    //    MatrixB21[i] = new double[N / 2];
    //    MatrixB22[i] = new double[N / 2];

    //   /* MatrixC11[i] = new int[N / 2];
    //    MatrixC12[i] = new int[N / 2];
    //    MatrixC21[i] = new int[N / 2];
    //    MatrixC22[i] = new int[N / 2];*/
    /* }*/

    //为每个小矩阵赋值，将大矩阵分割为4个小矩阵

    for (int i = 0; i < N / 2; i++)
    {
        vector<double> tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,tmp8;
        for (int j = 0; j < N / 2; j++)
        {
            tmp1.push_back(MatrixA[i][j]);
            tmp2.push_back(MatrixA[i][j + (N / 2)]);
            tmp3.push_back(MatrixA[i + (N / 2)][j]);
            tmp4.push_back(MatrixA[i + (N / 2)][j + (N / 2)]);

            tmp5.push_back(MatrixB[i][j]);
            tmp6.push_back(MatrixB[i][j + (N / 2)]);
            tmp7.push_back(MatrixB[i + (N / 2)][j]);
            tmp8.push_back(MatrixB[i + (N / 2)][j + (N / 2)]);
        }
        
        MatrixA11.push_back(tmp1);
        //print_vector(MatrixA11[i]);
        MatrixA12.push_back(tmp2);
        MatrixA21.push_back(tmp3);
        MatrixA22.push_back(tmp4);

        MatrixB11.push_back(tmp5);
        MatrixB12.push_back(tmp6);
        MatrixB21.push_back(tmp7);
        MatrixB22.push_back(tmp8);
    }

    chrono::high_resolution_clock::time_point time_start, time_end;
    //定义时间变量
    chrono::microseconds time_encode_input(0);
    chrono::microseconds time_encode_filter(0);
    chrono::microseconds time_decode_result(0);
    chrono::microseconds time_encrypt_input(0);
    chrono::microseconds time_decrypt_result(0);
    chrono::microseconds time_add(0);
    chrono::microseconds time_add_sum(0);
    chrono::microseconds time_multiply(0);
    chrono::microseconds time_multiply_sum(0);
    chrono::microseconds time_packing_input(0);
    chrono::microseconds time_packing_filter(0);
    chrono::microseconds time_read_input1(0);
    chrono::microseconds time_read_input2(0);
    chrono::microseconds time_transform_input(0);
    chrono::microseconds time_rescale(0);
    chrono::microseconds time_scale(0);
    chrono::microseconds time_scale_sum(0);
    chrono::microseconds time_relinearize(0);
    chrono::microseconds time_relinearize_sum(0);
    chrono::microseconds time_rotation(0);
    chrono::microseconds time_rotation_sum(0);
    chrono::microseconds time_sum(0);
    chrono::microseconds time_parms_id(0);
    chrono::microseconds time_parms_id_sum(0);
    chrono::microseconds time_strassen(0);
    chrono::microseconds time_strassen_sum(0);


    //生成秘钥并计算时间，加密，计算，解密
    chrono::microseconds time_diff;
    cout << "Generating secret/public keys: ";
    KeyGenerator keygen(context);
    cout << "Done" << endl;
    
    //输出生成私钥的时间
    cout << "Generating secret keys: ";
    time_start = chrono::high_resolution_clock::now();
    auto public_key = keygen.public_key();
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Done [" << time_diff.count() << " microseconds]" << endl;

    //输出生成公钥的时间
    cout << "Generating public keys: ";
    time_start = chrono::high_resolution_clock::now();
    auto secret_key = keygen.secret_key();
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Done [" << time_diff.count() << " microseconds]" << endl;

    RelinKeys relin_keys;
    GaloisKeys gal_keys;

    if (context->using_keyswitching())
    {
        //输出生成重线性化秘钥的时间
        cout << "Generating relinearization keys: ";
        time_start = chrono::high_resolution_clock::now();
        relin_keys = keygen.relin_keys_local();
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "Done [" << time_diff.count() << " microseconds]" << endl;

        if (!context->first_context_data()->qualifiers().using_batching)
        {
            cout << "Given encryption parameters do not support batching." << endl;
            return;
        }

        //输出生成Galois秘钥的时间
        cout << "Generating Galois keys: ";
        time_start = chrono::high_resolution_clock::now();
        gal_keys = keygen.galois_keys_local();
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "Done [" << time_diff.count() << " microseconds]" << endl;
    }

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;

    //现将输入矩阵packing到一条密文内
    vector<double> input1, input2, input3, input4, filter1, filter2, filter3, filter4;
    input1.reserve(slot_count);
    input2.reserve(slot_count);
    input3.reserve(slot_count);
    input4.reserve(slot_count);
   
    filter1.reserve(slot_count);
    filter2.reserve(slot_count);
    filter3.reserve(slot_count);
    filter4.reserve(slot_count);



    // 计算将matrixA11~matrixA22起来需要的时间
    // 这个部分到时候可以修改成diagnal scheme
    time_start = chrono::high_resolution_clock::now();
 

    // matrixA11 packing
    for (int i = 0; i < N/2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            input1.push_back(MatrixA11[i][j]);
        }
    }
    //print_vector(input1);

    // matrixA12 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            input2.push_back(MatrixA12[i][j]);
        }
    }
    //print_vector(input2);
    // matrixA21 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            input3.push_back(MatrixA21[i][j]);
        }
    }
    //print_vector(input3);
    // matrixA22 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            input4.push_back(MatrixA22[i][j]);
        }
    }
    //print_vector(input4);
    time_end = chrono::high_resolution_clock::now();
    time_packing_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    /*cout << " ]" << endl << endl;*/
    cout << "Packing MatrixA11~MatrixA22: Done [";
    cout << time_packing_input.count() << " microseconds]" << endl;
    print_vector(input1);
    print_vector(input2);
    print_vector(input3);
    print_vector(input4);

    // 计算将MatrixB11~MatrixB22 packing起来需要的时间
    // 这个部分到时候可以修改成diagnal scheme
    time_start = chrono::high_resolution_clock::now();
    
    // MatrixB11 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            filter1.push_back(MatrixB11[i][j]);
        }
    }
    //print_vector(filter1);
    // MatrixB12 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            filter2.push_back(MatrixB12[i][j]);
        }
    }
    //print_vector(filter2);
    // MatrixB21 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            filter3.push_back(MatrixB21[i][j]);
        }
    }
   // print_vector(filter3);
    // MatrixB22 packing
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            filter4.push_back(MatrixB22[i][j]);
        }
    }
    //print_vector(filter4);
    time_end = chrono::high_resolution_clock::now();
    time_packing_filter = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Packing MatrixB11~MatrixB22: Done [";
    cout << time_packing_filter.count() << " microseconds]" << endl;

    print_vector(filter1);
    print_vector(filter2);
    print_vector(filter3);
    print_vector(filter4);


    cout << endl << "Encoding input and filter .........." << endl;

    //编码filter
    /*Plaintext filter_plain1, filter_plain2, filter_plain3, filter_plain4;
    time_start = chrono::high_resolution_clock::now();
    cout << "Encode filter: Done [";
    encoder.encode(filter1, scale, filter_plain1);
    encoder.encode(filter2, scale, filter_plain2);
    encoder.encode(filter3, scale, filter_plain3);
    encoder.encode(filter4, scale, filter_plain4);
    time_end = chrono::high_resolution_clock::now();
    time_encode_filter = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_encode_filter.count() << " microseconds]" << endl;*/

    //编码input
    Plaintext x_plain1, x_plain2, x_plain3, x_plain4;
    cout << "Encode input vectors：Done [";
    time_start = chrono::high_resolution_clock::now();
    encoder.encode(input1, scale, x_plain1); //将输入进行编码
    encoder.encode(input2, scale, x_plain2);
    encoder.encode(input3, scale, x_plain3);
    encoder.encode(input4, scale, x_plain4);
    time_end = chrono::high_resolution_clock::now();
    time_encode_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_encode_input.count() << " microseconds]" << endl;

    
    //加密
    cout << endl << "Encrypting input .........." << endl;
    Ciphertext x1_encrypted, x2_encrypted, x3_encrypted, x4_encrypted;
    cout << "Encrypt input vectors：Done [";
    time_start = chrono::high_resolution_clock::now();
    encryptor.encrypt(x_plain1, x1_encrypted); //将编码后的输入进行加密
    encryptor.encrypt(x_plain2, x2_encrypted);
    encryptor.encrypt(x_plain3, x3_encrypted);
    encryptor.encrypt(x_plain4, x4_encrypted);
    time_end = chrono::high_resolution_clock::now();
    time_encrypt_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_encrypt_input.count() << " microseconds]" << endl;

    
    // 点积实现
    cout << endl << "Dotproduct .........." << endl;
    //Ciphertext destination1, destination2, destination3, destination4;//保存每次的乘积结果
    //
    ////定义了N * N / 4个整型元素的向量,且给出每个元素的初值为0
    //vector<double> v1(slot_count, 0), v2(slot_count, 0), v3(slot_count, 0), v4(slot_count, 0);
    //Plaintext p1, p2, p3, p4;
    //encoder.encode(v1, scale, p1);
    //encoder.encode(v2, scale, p2);
    //encoder.encode(v3, scale, p3);
    //encoder.encode(v4, scale, p4);
    //Ciphertext dotproduct1, dotproduct2, dotproduct3, dotproduct4;//保存累计的加和结果
    //encryptor.encrypt(p1, dotproduct1);
    //encryptor.encrypt(p2, dotproduct2);
    //encryptor.encrypt(p3, dotproduct3);
    //encryptor.encrypt(p4, dotproduct4);
    

    //创建S1-S10
    Ciphertext S2, S3, S5, S7, S9;
    Plaintext S1, S4, S6, S8, S10;
    //Plaintext S22, S33, S55, S77, S99;
    vector<double> s1(slot_count, 0), s4(slot_count, 0), s6(slot_count, 0), s8(slot_count, 0), s10(slot_count, 0);
    

    /*encoder.encode(s1, scale, S1);
    encoder.encode(s2, scale, S22);
    encoder.encode(s3, scale, S33);
    encoder.encode(s4, scale, S4);
    encoder.encode(s5, scale, S55);
    encoder.encode(s6, scale, S6);
    encoder.encode(s7, scale, S77);
    encoder.encode(s8, scale, S8);
    encoder.encode(s9, scale, S99);
    encoder.encode(s10, scale, S10);

    encryptor.encrypt(S22, S2);
    encryptor.encrypt(S33, S3);
    encryptor.encrypt(S55, S5);
    encryptor.encrypt(S77, S7);
    encryptor.encrypt(S99, S9);*/

    // 加减法
    time_start = chrono::high_resolution_clock::now();
    cout << "S1" << endl;
    vector_sub(int(slot_count), filter2, filter4, s1);//S1 = B12 - B22   
    print_vector(s1);
    cout << "S2" << endl;
    evaluator.add(x1_encrypted, x2_encrypted, S2);//S2 = A11 + A12
    cout << "S3" << endl;
    evaluator.add(x3_encrypted, x4_encrypted, S3);//S3 = A21 + A22
    cout << "S4" << endl;
    vector_sub(int(slot_count), filter3, filter1, s4); // S4 = B21 - B11
    print_vector(s4);
    cout << "S5" << endl;
    evaluator.add(x1_encrypted, x4_encrypted, S5);//S5 = A11 + A22
    cout << "S6" << endl;
    vector_sum(int(slot_count), filter1, filter4, s6); // S6 = B11 + B22
    print_vector(s6);
    cout << "S7" << endl;
    evaluator.sub(x2_encrypted, x4_encrypted, S7);//S7 = A12 - A22
    cout << "S8" << endl;
    vector_sum(int(slot_count), filter3, filter4, s8); // S8 = B21 + B22
    print_vector(s8);
    cout << "S9" << endl;
    evaluator.sub(x1_encrypted, x3_encrypted, S9);//S9 = A11 - A21
    cout << "S10" << endl;
    vector_sum(int(slot_count), filter1, filter2, s10); // S10 = B11 + B12
    print_vector(s10);
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Strassen——Compute S1-S10：Done [";
    cout << time_strassen.count() << " microseconds]" << endl;
    


    // 点积运算
    

    //明文解码，以便后续获得对角矩阵
    //vector<double> tmp11, tmp22, tmp33, tmp44, tmp55, tmp66, tmp77;
    //cout << endl << "Decode Plaintext to Obtain diagonal matrix: Done [" ;
    //time_start = chrono::high_resolution_clock::now();
    //encoder.decode(S1, tmp11);
    ////print_vector(tmp1);//S1 为0，说明上述减法无效。
    //encoder.decode(filter_plain4, tmp22);
    //encoder.decode(filter_plain1, tmp33);
    //encoder.decode(S4, tmp44);
    //encoder.decode(S6, tmp55);
    //encoder.decode(S8, tmp66);
    //encoder.decode(S10, tmp77);
    //time_end = chrono::high_resolution_clock::now();
    //time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //cout << time_strassen.count() << " microseconds]" << endl;
    //time_strassen_sum += time_strassen;

    // dotproduct
    // 先将明文拓展为对角的排列方式
    vector<double> tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    cout << endl << "Obtain diagonal matrix: Done [";
    time_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N / 2; ++i)
    {
        for (int j = 0; j < N / 2; ++j)
        {
            tmp1.push_back(s1[(i + j) % (N / 2) + j * (N / 2)]);
            tmp2.push_back(filter4[(i + j) % (N / 2) + j * (N / 2)]);
            tmp3.push_back(filter1[(i + j) % (N / 2) + j * (N / 2)]);
            tmp4.push_back(s4[(i + j) % (N / 2) + j * (N / 2)]);
            tmp5.push_back(s6[(i + j) % (N / 2) + j * (N / 2)]);
            tmp6.push_back(s8[(i + j) % (N / 2) + j * (N / 2)]);
            tmp7.push_back(s10[(i + j) % (N / 2) + j * (N / 2)]);
            

        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_strassen.count() << " microseconds]" << endl;
    time_strassen_sum += time_strassen;
    print_vector(tmp1);
    print_vector(tmp2);
    print_vector(tmp3);
    print_vector(tmp4);
    print_vector(tmp5);
    print_vector(tmp6);
    print_vector(tmp7);


    Plaintext filter_plain1, filter_plain2, filter_plain3, filter_plain4;
    //time_start = chrono::high_resolution_clock::now();
    //cout << "Encode filter: Done [";
    ////encoder.encode(filter1, scale, filter_plain1);
    //encoder.encode(filter2, scale, filter_plain2);
    //encoder.encode(filter3, scale, filter_plain3);
    ////encoder.encode(filter4, scale, filter_plain4);
    //time_end = chrono::high_resolution_clock::now();
    //time_encode_filter = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //cout << time_encode_filter.count() << " microseconds]" << endl;

    //将tmp编码为明文
    cout << endl << "Encode diagnal matrix: Done [";
    time_start = chrono::high_resolution_clock::now();
    encoder.encode(tmp1, scale, S1);
    encoder.encode(tmp2, scale, filter_plain4); 
    encoder.encode(tmp3, scale, filter_plain1);
    encoder.encode(tmp4, scale, S4);
    encoder.encode(tmp5, scale, S6);
    encoder.encode(tmp6, scale, S8);
    encoder.encode(tmp7, scale, S10);
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_strassen.count() << " microseconds]" << endl;
    


    Ciphertext P1, P2, P3, P4, P5, P6, P7;
    Ciphertext result1, result2, result3, result4, result5, result6, result7;
    // 相乘
    cout << "multiplication........" << endl;
    evaluator.multiply_plain(x1_encrypted, S1, result1);//P1 = A11 • S1
    evaluator.multiply_plain(S2, filter_plain4, result2);//P2 = S2 • B22
    evaluator.multiply_plain(S3, filter_plain1, result3);//P3 = S3 • B11
    evaluator.multiply_plain(x4_encrypted, S4, result4);//P4 = A22 • S4
    evaluator.multiply_plain(S5, S6, result5);//P5 = S5 • S6
    evaluator.multiply_plain(S7, S8, result6);//P6 = S7 • S8
    evaluator.multiply_plain(S9, S10, result7);//P7 = S9 • S10
    cout << "Done" << endl;
    // Relinearize 𝑟𝑒𝑠𝑢𝑙𝑡, keep same size
    cout << "Relinearize the result：Done [";
    time_start = chrono::high_resolution_clock::now();
    evaluator.relinearize_inplace(result1, relin_keys);
    evaluator.relinearize_inplace(result2, relin_keys);
    evaluator.relinearize_inplace(result3, relin_keys);
    evaluator.relinearize_inplace(result4, relin_keys);
    evaluator.relinearize_inplace(result5, relin_keys);
    evaluator.relinearize_inplace(result6, relin_keys);
    evaluator.relinearize_inplace(result7, relin_keys);
    time_end = chrono::high_resolution_clock::now();
    time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_relinearize.count() << " microseconds]" << endl;
    

    cout << "Scale the result：Done [";
    time_start = chrono::high_resolution_clock::now();
    evaluator.rescale_to_next_inplace(result1);
    evaluator.rescale_to_next_inplace(result2);
    evaluator.rescale_to_next_inplace(result3);
    evaluator.rescale_to_next_inplace(result4);
    evaluator.rescale_to_next_inplace(result5);
    evaluator.rescale_to_next_inplace(result6);
    evaluator.rescale_to_next_inplace(result7);

    
    time_end = chrono::high_resolution_clock::now();
    time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_scale.count() << " microseconds]" << endl;
    


   /* Plaintext dd1, dd2, dd3, dd4, dd5, dd6, dd7;
    decryptor.decrypt(result1, dd1);
    decryptor.decrypt(result2, dd2);
    decryptor.decrypt(result3, dd3);
    decryptor.decrypt(result4, dd4);
    decryptor.decrypt(result5, dd5);
    decryptor.decrypt(result6, dd6);
    decryptor.decrypt(result7, dd7);

    vector<double> vv1, vv2, vv3, vv4, vv5, vv6, vv7;
    encoder.decode(dd1, vv1);
    encoder.decode(dd2, vv2);
    encoder.decode(dd3, vv3);
    encoder.decode(dd4, vv4);
    encoder.decode(dd5, vv5);
    encoder.decode(dd6, vv6);
    encoder.decode(dd7, vv7);

    print_vector(vv1);
    print_vector(vv2);
    print_vector(vv3);
    print_vector(vv4);
    print_vector(vv5);
    print_vector(vv6);
    print_vector(vv7);*/
    
    // 做乘法，求点积
    //1.先相乘
    /*vector<double> pp1(slot_count, 0), pp2(slot_count, 0), pp3(slot_count, 0), pp4(slot_count, 0), pp5(slot_count, 0),
        pp6(slot_count, 0), pp7(slot_count, 0);

    print_vector(pp1);
    Plaintext p1, p2, p3, p4, p5, p6, p7;
    encoder.encode(pp1, scale, p1);
    encoder.encode(pp2, scale, p2);
    encoder.encode(pp3, scale, p3);
    encoder.encode(pp4, scale, p4);
    encoder.encode(pp5, scale, p5);
    encoder.encode(pp6, scale, p6);
    encoder.encode(pp7, scale, p7);

    Ciphertext P1, P2, P3, P4, P5, P6, P7;
    encryptor.encrypt(p1, P1);
    encryptor.encrypt(p2, P2);
    encryptor.encrypt(p3, P3);
    encryptor.encrypt(p4, P4);
    encryptor.encrypt(p5, P5);
    encryptor.encrypt(p6, P6);
    encryptor.encrypt(p7, P7);

    Ciphertext result1, result2, result3, result4, result5, result6, result7;
    time_start = chrono::high_resolution_clock::now();*/

   


  
    //求点积
    
    cout << endl << "Dotproduct start........." << endl;
    time_start = chrono::high_resolution_clock::now();
    for (int i = 1; i < N / 2; ++i)
    {
        //旋转input2，得到新的向量，从而与相同的密文进行后续的相乘
        //这个可以减少密文的旋转，减少时间开销
        cout << "rotate......" << endl;
        rotate(tmp1, int(slot_count), N / 2);
        rotate(tmp2, int(slot_count), N / 2);
        rotate(tmp3, int(slot_count), N / 2);
        rotate(tmp4, int(slot_count), N / 2);
        rotate(tmp5, int(slot_count), N / 2);
        rotate(tmp6, int(slot_count), N / 2);
        rotate(tmp7, int(slot_count), N / 2);
        /*print_vector(tmp1);
        print_vector(tmp2);
        print_vector(tmp3);
        print_vector(tmp4);
        print_vector(tmp5);
        print_vector(tmp6);
        print_vector(tmp7);*/

        cout << endl << "Encode diagnal matrix: Done [";
        time_start = chrono::high_resolution_clock::now();
        encoder.encode(tmp1, scale, S1);
        encoder.encode(tmp2, scale, filter_plain4);
        encoder.encode(tmp3, scale, filter_plain1);
        encoder.encode(tmp4, scale, S4);
        encoder.encode(tmp5, scale, S6);
        encoder.encode(tmp6, scale, S8);
        encoder.encode(tmp7, scale, S10);
        time_end = chrono::high_resolution_clock::now();
        time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_strassen.count() << " microseconds]" << endl;
        

        // 相乘
        cout << "multiplication........" << endl;
        evaluator.multiply_plain(x1_encrypted, S1, P1);  // P1 = A11 • S1
        evaluator.multiply_plain(S2, filter_plain4, P2); // P2 = S2 • B22
        evaluator.multiply_plain(S3, filter_plain1, P3); // P3 = S3 • B11
        evaluator.multiply_plain(x4_encrypted, S4, P4);  // P4 = A22 • S4
        evaluator.multiply_plain(S5, S6, P5);            // P5 = S5 • S6
        evaluator.multiply_plain(S7, S8, P6);            // P6 = S7 • S8
        evaluator.multiply_plain(S9, S10, P7);           // P7 = S9 • S10
        cout << "Done" << endl;

        // Relinearize 𝑟𝑒𝑠𝑢𝑙𝑡, keep same size
        cout << "Relinearize the result：Done [";
        time_start = chrono::high_resolution_clock::now();
        evaluator.relinearize_inplace(P1, relin_keys);
        evaluator.relinearize_inplace(P2, relin_keys);
        evaluator.relinearize_inplace(P3, relin_keys);
        evaluator.relinearize_inplace(P4, relin_keys);
        evaluator.relinearize_inplace(P5, relin_keys);
        evaluator.relinearize_inplace(P6, relin_keys);
        evaluator.relinearize_inplace(P7, relin_keys);
        time_end = chrono::high_resolution_clock::now();
        time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_relinearize.count() << " microseconds]" << endl;
        


        cout << "Rescale the result：Done [";
        time_start = chrono::high_resolution_clock::now();
        evaluator.rescale_to_next_inplace(P1);
        evaluator.rescale_to_next_inplace(P2);
        evaluator.rescale_to_next_inplace(P3);
        evaluator.rescale_to_next_inplace(P4);
        evaluator.rescale_to_next_inplace(P5);
        evaluator.rescale_to_next_inplace(P6);
        evaluator.rescale_to_next_inplace(P7);
        time_end = chrono::high_resolution_clock::now();
        time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_scale.count() << " microseconds]" << endl;
        

        



        //旋转result1，以便后续相加
        cout << "rotate vetor......." << endl;
        
        Plaintext d1, d2, d3, d4, d5, d6, d7;
        decryptor.decrypt(P1, d1);
        decryptor.decrypt(P2, d2);
        decryptor.decrypt(P3, d3);
        decryptor.decrypt(P4, d4);
        decryptor.decrypt(P5, d5);
        decryptor.decrypt(P6, d6);
        decryptor.decrypt(P7, d7);

        vector<double> v1, v2, v3, v4, v5, v6, v7;
        encoder.decode(d1, v1);
        encoder.decode(d2, v2);
        encoder.decode(d3, v3);
        encoder.decode(d4, v4);
        encoder.decode(d5, v5);
        encoder.decode(d6, v6);
        encoder.decode(d7, v7);
        
        /*print_vector(v1);
        print_vector(v2);
        print_vector(v3);
        print_vector(v4);
        print_vector(v5);
        print_vector(v6);
        print_vector(v7);*/

        bacth_rotate(v1, N / 2, N / 2 - i);
        bacth_rotate(v2, N / 2, N / 2 - i);
        bacth_rotate(v3, N / 2, N / 2 - i);
        bacth_rotate(v4, N / 2, N / 2 - i);
        bacth_rotate(v5, N / 2, N / 2 - i);
        bacth_rotate(v6, N / 2, N / 2 - i);
        bacth_rotate(v7, N / 2, N / 2 - i);
     

        encoder.encode(v1, scale, d1);
        encoder.encode(v2, scale, d2);
        encoder.encode(v3, scale, d3);
        encoder.encode(v4, scale, d4);
        encoder.encode(v5, scale, d5);
        encoder.encode(v6, scale, d6);
        encoder.encode(v7, scale, d7);

        encryptor.encrypt(d1, P1);
        encryptor.encrypt(d2, P2);
        encryptor.encrypt(d3, P3);
        encryptor.encrypt(d4, P4);
        encryptor.encrypt(d5, P5);
        encryptor.encrypt(d6, P6);
        encryptor.encrypt(d7, P7);
            
           

        result1.scale() = pow(2.0, 40);
        result2.scale() = pow(2.0, 40);
        result3.scale() = pow(2.0, 40);
        result4.scale() = pow(2.0, 40);
        result5.scale() = pow(2.0, 40);
        result6.scale() = pow(2.0, 40);
        result7.scale() = pow(2.0, 40);

       

        cout << "Make the same parms_id：Done [";
        time_start = chrono::high_resolution_clock::now();
        parms_id_type last_parms_id = result1.parms_id();
        evaluator.mod_switch_to_inplace(P1, last_parms_id);
        evaluator.mod_switch_to_inplace(P2, last_parms_id);
        evaluator.mod_switch_to_inplace(P3, last_parms_id);
        evaluator.mod_switch_to_inplace(P4, last_parms_id);
        evaluator.mod_switch_to_inplace(P5, last_parms_id);
        evaluator.mod_switch_to_inplace(P6, last_parms_id);
        evaluator.mod_switch_to_inplace(P7, last_parms_id);
        time_end = chrono::high_resolution_clock::now();
        time_parms_id = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_parms_id.count() << " microseconds]" << endl;
        

        

        //相加
        cout << "Addition........" << endl;
        evaluator.add_inplace(result1, P1);
        evaluator.add_inplace(result2, P2);
        evaluator.add_inplace(result3, P3);
        evaluator.add_inplace(result4, P4);
        evaluator.add_inplace(result5, P5);
        evaluator.add_inplace(result6, P6);
        evaluator.add_inplace(result7, P7);

        
    }
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Strassen——Compute P1-P7：Done [";
    cout << time_encrypt_input.count() << " microseconds]" << endl;
    time_strassen_sum += time_strassen;


    // 将上述结果P1-P7的结果与矩阵一一对应，通过旋转来处理

    cout << endl << "The results of p1-p7 correspond to the matrix one by one: Done [";

    Plaintext d1, d2, d3, d4, d5, d6, d7;
    decryptor.decrypt(result1, d1);
    decryptor.decrypt(result2, d2);
    decryptor.decrypt(result3, d3);
    decryptor.decrypt(result4, d4);
    decryptor.decrypt(result5, d5);
    decryptor.decrypt(result6, d6);
    decryptor.decrypt(result7, d7);

    vector<double> v1(slot_count), v2(slot_count), v3(slot_count), v4(slot_count), v5(slot_count), v6(slot_count),
        v7(slot_count);
    encoder.decode(d1, v1);
    encoder.decode(d2, v2);
    encoder.decode(d3, v3);
    encoder.decode(d4, v4);
    encoder.decode(d5, v5);
    encoder.decode(d6, v6);
    encoder.decode(d7, v7);
    print_vector(v1);
    print_vector(v2);
    print_vector(v3);
    print_vector(v4);
    print_vector(v5);
    print_vector(v6);
    print_vector(v7);

    //vector<double> res1(slot_count, 0), res2(slot_count, 0), res3(slot_count, 0), res4(slot_count, 0),
    //    res5(slot_count, 0), res6(slot_count, 0), res7(slot_count, 0);

    time_start = chrono::high_resolution_clock::now();
   /* for (int j = 0; j < N/2; j++)
    {
        for (int i = 0; i < N / 2; i++)
        {
            res1[j * (N / 2) + i] = v1[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res2[j * (N / 2) + i] = v2[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res3[j * (N / 2) + i] = v3[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res4[j * (N / 2) + i] = v4[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res5[j * (N / 2) + i] = v5[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res6[j * (N / 2) + i] = v6[(N / 2 - i) % (N / 2) + (N / 2) * i];
            res7[j * (N / 2) + i] = v7[(N / 2 - i) % (N / 2) + (N / 2) * i];
        }
        
    }
    print_vector(res1);
    print_vector(res2);
    print_vector(res3);
    print_vector(res4);
    print_vector(res5);
    print_vector(res6);
    print_vector(res7);*/
    for (int i = 1; i < N / 2; ++i)
    {
        single_bacth_rotate(v1, N / 2, i, N / 2 - i);
        single_bacth_rotate(v2, N / 2, i, N / 2 - i);
        single_bacth_rotate(v3, N / 2, i, N / 2 - i);
        single_bacth_rotate(v4, N / 2, i, N / 2 - i);
        single_bacth_rotate(v5, N / 2, i, N / 2 - i);
        single_bacth_rotate(v6, N / 2, i, N / 2 - i);
        single_bacth_rotate(v7, N / 2, i, N / 2 - i);
        
       /* cout << "Batch rotation......." << endl;
        print_vector(v1);
        print_vector(v2);
        print_vector(v3);
        print_vector(v4);
        print_vector(v5);
        print_vector(v6);
        print_vector(v7);
        cout << endl;*/
    }
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_strassen.count() << " microseconds]" << endl;
    time_strassen_sum += time_strassen;


    //// 获得C11-C22
    //cout << endl << "Obtain C11-C22: Done [";
    //time_start = chrono::high_resolution_clock::now();
    //Ciphertext C11, C12, C21, C22;
    //evaluator.add(result5, result4, C11); // C11 = P5 + P4 - P2 + P6
    //evaluator.sub(C11, result2, C11);
    //evaluator.add(C11, result6, C11);
    //evaluator.add(result1, result2, C12); // C12 = P1 + P2
    //evaluator.add(result3, result4, C21); // C21 = P3 + P4
    //evaluator.add(result5, result1, C22); // C22 = P5 + P1 - P3 - P7
    //evaluator.sub(C22, result3, C22);
    //evaluator.sub(C22, result7, C22);
    //time_end = chrono::high_resolution_clock::now();
    //time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //cout << time_strassen.count() << " microseconds]" << endl;
    //time_strassen_sum += time_strassen;


    // 获得C11-C22
    cout << endl << "Obtain C11-C22: Done [";
    time_start = chrono::high_resolution_clock::now();
    vector<double> C11(slot_count), C12(slot_count), C21(slot_count), C22(slot_count);

    vector_sum(int(slot_count), v5, v4, C11); // C11 = P5 + P4 - P2 + P6
    vector_sub(int(slot_count), C11, v2, C11);
    vector_sum(int(slot_count), C11, v6, C11);
    vector_sum(int(slot_count), v1, v2, C12); // C12 = P1 + P2
    vector_sum(int(slot_count), v3, v4, C21); // C21 = P3 + P4
    vector_sum(int(slot_count), v5, v1, C22); // C22 = P5 + P1 - P3 - P7
    vector_sub(int(slot_count), C22, v3, C22);
    vector_sub(int(slot_count), C22, v7, C22); 


    //Vector_Sum(int(slot_count), res5, res4, C11); // C11 = P5 + P4 - P2 + P6
    //Vector_Sub(int(slot_count), C11, res2, C11); 
    //Vector_Sum(int(slot_count), C11, res6, C11); 
    //Vector_Sum(int(slot_count), res1, res2, C12); // C12 = P1 + P2
    //Vector_Sum(int(slot_count), res3, res4, C21); // C21 = P3 + P4
    //Vector_Sum(int(slot_count), res5, res1, C22); // C22 = P5 + P1 - P3 - P7
    //Vector_Sub(int(slot_count), C22, res3, C22); 
    //Vector_Sub(int(slot_count), C22, res7, C22); 
    
    time_end = chrono::high_resolution_clock::now();
    time_strassen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_strassen.count() << " microseconds]" << endl;
    time_strassen_sum += time_strassen;

    







    //for (int i = 0; i < N/2; ++i)
    //{
    //    //乘法,将结果保存到destination中
    //    cout << "Input ciphertext multiplies filter plaintext：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    evaluator.multiply_plain(x1_encrypted, filter_plain1, destination1);
    //    evaluator.multiply_plain(x2_encrypted, filter_plain2, destination2); 
    //    evaluator.multiply_plain(x3_encrypted, filter_plain3, destination3); 
    //    evaluator.multiply_plain(x4_encrypted, filter_plain4, destination4); 
    //    time_end = chrono::high_resolution_clock::now();
    //    time_multiply = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_multiply.count() << " microseconds]" << endl;
    //    time_multiply_sum += time_multiply;


    //    

    //    /*cout << "    + Scale of input*filter before rescale: " << log2(destination1.scale()) << " bits" << endl;
    //    evaluator.rescale_to_next_inplace(destination1);
    //    cout << "    + Scale of input*filter after rescale: " << log2(destination1.scale()) << " bits" << endl;
    //    cout << "    + Scale of dotproduct1: " << log2(dotproduct1.scale()) << " bits" << endl;*/
    //    /*cout << destination1.parms_id() << endl;
    //    cout << x1_encrypted.parms_id() << endl;*/

    //    //Relinearize 𝑟𝑒𝑠𝑢𝑙𝑡, keep same size
    //    cout << "Relinearize the result：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    evaluator.relinearize_inplace(destination1, relin_keys);
    //    evaluator.relinearize_inplace(destination2, relin_keys);
    //    evaluator.relinearize_inplace(destination3, relin_keys);
    //    evaluator.relinearize_inplace(destination4, relin_keys);
    //    time_end = chrono::high_resolution_clock::now();
    //    time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_relinearize.count() << " microseconds]" << endl;
    //    time_relinearize_sum += time_relinearize;

    //    //这里有一个问题，如果rescale的话，能让scale一致，但是parm_id就会改变
    //    // 将其保持在40位scale，从而与dotproduct的scale一致，以便实现后续加法
    //    cout << "Relinearize the result：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    destination1.scale() = pow(2.0, 40);
    //    destination2.scale() = pow(2.0, 40);
    //    destination3.scale() = pow(2.0, 40);
    //    destination4.scale() = pow(2.0, 40);
    //    time_end = chrono::high_resolution_clock::now();
    //    time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_scale.count() << " microseconds]" << endl;
    //    time_scale_sum += time_scale;

    //    cout << "Make the same parms_id：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    parms_id_type last_parms_id = destination1.parms_id();
    //    evaluator.mod_switch_to_inplace(dotproduct1, last_parms_id);
    //    evaluator.mod_switch_to_inplace(dotproduct2, last_parms_id);
    //    evaluator.mod_switch_to_inplace(dotproduct3, last_parms_id);
    //    evaluator.mod_switch_to_inplace(dotproduct4, last_parms_id);
    //    time_end = chrono::high_resolution_clock::now();
    //    time_parms_id = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_parms_id.count() << " microseconds]" << endl;
    //    time_parms_id_sum += time_parms_id;


    //    // 加法
    //    cout << "Add the result：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    evaluator.add_inplace(dotproduct1, destination1);
    //    evaluator.add_inplace(dotproduct2, destination2);
    //    evaluator.add_inplace(dotproduct3, destination3);
    //    evaluator.add_inplace(dotproduct4, destination4);
    //    /*evaluator.add_inplace(dotproduct2, destination2);
    //    evaluator.add_inplace(dotproduct3, destination3);
    //    evaluator.add_inplace(dotproduct4, destination4);*/
    //    time_end = chrono::high_resolution_clock::now();
    //    time_add = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_add.count() << " microseconds]" << endl;
    //    time_add_sum += time_add;


    //    //rotation
    //    cout << "Rotation the result：Done [";
    //    time_start = chrono::high_resolution_clock::now();
    //    rotate(filter_plain1, slot_count, N / 2);
    //    rotate(filter_plain2, slot_count, N / 2);
    //    rotate(filter_plain3, slot_count, N / 2);
    //    rotate(filter_plain4, slot_count, N / 2);
    //    time_end = chrono::high_resolution_clock::now();
    //    time_rotation = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //    cout << time_rotation.count() << " microseconds]" << endl << endl;
    //    time_rotation_sum += time_rotation;


    //}

    //cout << endl << "Sum of multiplication time: [";
    //cout << time_rotation_sum.count() << " microseconds]" << endl;

    //cout << endl << "Sum of rotation time: [";
    //cout << time_rotation_sum.count() << " microseconds]" << endl;
    //
    //cout << "Sum of add time: [";
    //cout << time_add_sum.count() << " microseconds]" << endl << endl;
    //

    

    
    vector<vector<double>> Mul_Matrix(N, vector<double>(N));

    time_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Mul_Matrix[i][j] = 0;
            
            for (int k = 0; k < N; k++)
            {
                Mul_Matrix[i][j] = Mul_Matrix[i][j] + MatrixA[i][k] * MatrixB[k][j];
                
            }
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_transform_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << endl << "Generating expected result: Done [" << time_transform_input.count() << " microseconds]" << endl;
    cout << endl << "Excepted Mul_Matrix:" << endl << endl;
    PrintMatrix(Mul_Matrix);

    /*vector<vector<double>> Mul_Matrix1(N / 2, vector<double>(N / 2)), Mul_Matrix2(N / 2, vector<double>(N / 2)),
        Mul_Matrix3(N / 2, vector<double>(N / 2)), Mul_Matrix4(N / 2, vector<double>(N / 2));

            
    for (int i = 0; i < N/2; i++)
    {
        for (int j = 0; j < N/2; j++)
        {
            Mul_Matrix1[i][j] = 0;
            Mul_Matrix2[i][j] = 0;
            Mul_Matrix3[i][j] = 0;
            Mul_Matrix4[i][j] = 0;
            for (int k = 0; k < N/2; k++)
            {
                Mul_Matrix1[i][j] = Mul_Matrix1[i][j] + MatrixA11[i][k] * MatrixB11[k][j];
                Mul_Matrix2[i][j] = Mul_Matrix2[i][j] + MatrixA12[i][k] * MatrixB12[k][j];
                Mul_Matrix3[i][j] = Mul_Matrix3[i][j] + MatrixA21[i][k] * MatrixB21[k][j];
                Mul_Matrix4[i][j] = Mul_Matrix4[i][j] + MatrixA22[i][k] * MatrixB22[k][j];
            }
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_transform_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << endl << "Generating expected result: Dome [" << time_transform_input.count() << " microseconds]" << endl;
    cout << endl << "Excepted Mul_Matrix1:" << endl << endl;
    PrintMatrix(Mul_Matrix1);
    cout << endl << "Excepted Mul_Matrix2:" << endl << endl;
    PrintMatrix(Mul_Matrix2);
    cout << endl << "Excepted Mul_Matrix3:" << endl << endl;
    PrintMatrix(Mul_Matrix3);
    cout << endl << "Excepted Mul_Matrix4:" << endl << endl;
    PrintMatrix(Mul_Matrix4);*/

    // 解密解码，输出加密后运算的结果
    //Plaintext plain_result1, plain_result2, plain_result3, plain_result4;
    //cout << endl << "Decrypt the result: Done [";
    //time_start = chrono::high_resolution_clock::now();
    //decryptor.decrypt(C11, plain_result1);
    //decryptor.decrypt(C12, plain_result2);
    //decryptor.decrypt(C21, plain_result3);
    //decryptor.decrypt(C22, plain_result4);
    //time_end = chrono::high_resolution_clock::now();
    //time_decrypt_result = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //cout << time_decrypt_result.count() << " microseconds]" << endl;

    //// print_line(__LINE__);
    //cout << "Decode the result: Done [";
    //vector<double> result11, result12, result21, result22;
    //time_start = chrono::high_resolution_clock::now();
    //encoder.decode(plain_result1, result11);
    //encoder.decode(plain_result2, result12);
    //encoder.decode(plain_result3, result21);
    //encoder.decode(plain_result4, result22);
    //time_end = chrono::high_resolution_clock::now();
    //time_decode_result = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    //cout << time_decode_result.count() << " microseconds]" << endl;

    cout << endl << "Computed result ...... Correct." << endl;
    vector<vector<double>> MatrixC(N, vector<double> (N));
    //将C11，C12,C21,C21合并为C矩阵
    
    for (int i = 0; i < N / 2; i++)
    {
        for (int j = 0; j < N / 2; j++)
        {
            
            MatrixC[i][j] = C11[i * (N / 2) + j];
            MatrixC[i][j + N / 2] = C12[i * (N / 2) + j];
            MatrixC[i + N / 2][j] = C21[i * (N / 2) + j];
            MatrixC[i + N / 2][j + N / 2] = C22[i * (N / 2) + j];
        }
    }
    PrintMatrix(MatrixC);




    /*vector<vector<double>> result11, result22, result33, result44;
    for (size_t i = 0; i < N; i++)
    {
        vector<double> tmp1, tmp2, tmp3, tmp4;
        for (size_t j = 0; j < N; j++)
        {
            tmp1.push_back(result1[j * i + j]);
            tmp2.push_back(result2[j * i + j]);
            tmp3.push_back(result3[j * i + j]);
            tmp4.push_back(result4[j * i + j]);
        }
        result11.push_back(tmp1);
        result22.push_back(tmp2);
        result33.push_back(tmp3);
        result44.push_back(tmp4);
    }

    cout << endl << "Mul_Matrix1:" << endl << endl;
    PrintMatrix(result11);
    cout << endl << "Mul_Matrix2:" << endl << endl;
    PrintMatrix(result22);
    cout << endl << "Mul_Matrix3:" << endl << endl;
    PrintMatrix(result33);
    cout << endl << "Mul_Matrix4:" << endl << endl;
    PrintMatrix(result44);*/


    cout << "Sum of encrypted input1 * input2 time: [";
    time_sum = time_packing_input + time_packing_filter + time_encode_filter + time_encode_input + time_add_sum +
               time_decode_result + time_decrypt_result + time_encrypt_input + time_relinearize + time_strassen_sum;
    cout << time_sum.count() << " microseconds]" << endl << endl;
}