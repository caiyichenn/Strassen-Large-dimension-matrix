// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "examples.h"

using namespace std;
using namespace seal;

void example_ckks_dotproduct()
{
    print_example_banner("Example: CKKS Basics");

    /*
    In this example we demonstrate evaluating a polynomial function
    在这个例子中，我们演示了如何在加密的浮点输入数据x上计算一个多项式函数PI*x^3 + 0.4*x + 1，
    在间隔为[0,1]的一组4096个等距点上。这个例子展示了CKKS方案的许多主要特性，但也展示了使用它时的挑战。

        PI*x^3 + 0.4*x + 1

    on encrypted floating-point input data x for a set of 4096 equidistant points
    in the interval [0, 1]. This example demonstrates many of the main features
    of the CKKS scheme, but also the challenges in using it.

    We start by setting up the CKKS scheme.
    */
    chrono::high_resolution_clock::time_point time_start, time_end;
    EncryptionParameters parms(scheme_type::CKKS);

    /*
    We saw in `2_encoders.cpp' that multiplication in CKKS causes scales
    in ciphertexts to grow. The scale of any ciphertext must not get too close
    to the total size of coeff_modulus, or else the ciphertext simply runs out of
    room to store the scaled-up plaintext. The CKKS scheme provides a `rescale'
    functionality that can reduce the scale, and stabilize the scale expansion.
    我们在' 2_encoders.cpp'中看到CKKS中的乘法会导致密文的规模增长。
    任何密文的规模都不能太接近coeff_modulus的总大小，否则密文就会耗尽存储放大后
    的明文的空间。CKKS方案提供了“缩放rescale”功能，可以减少规模，并稳定规模扩展。

    Rescaling is a kind of modulus switch operation (recall `3_levels.cpp').
    As modulus switching, it removes the last of the primes from coeff_modulus,
    but as a side-effect it scales down the ciphertext by the removed prime.
    Usually we want to have perfect control over how the scales are changed,
    which is why for the CKKS scheme it is more common to use carefully selected
    primes for the coeff_modulus.
    缩放是一种模量切换操作(回顾'3_levels.cpp')。作为模量切换，它从coeff_modulus中删除
    了最后一个素数，但作为一个副作用，它通过删除的素数来缩小密文。通常我们想要完美
    地控制尺度的变化，这就是为什么对于CKKS方案，更常见的是使用精心挑选的质数作为系数模。

    More precisely, suppose that the scale in a CKKS ciphertext is S, and the
    last prime in the current coeff_modulus (for the ciphertext) is P. Rescaling
    to the next level changes the scale to S/P, and removes the prime P from the
    coeff_modulus, as usual in modulus switching. The number of primes limits
    how many rescalings can be done, and thus limits the multiplicative depth of
    the computation.
    更准确地说，假设CKKS密文的尺度是S，而当前的系数模(对于密文)中的最后一个素数是P。
    缩放到下一层时，将尺度改变为S/P，并从系数模中删除素数P，就像模切换中通常的情况一样。
    质数的数量限制了缩放的次数，从而限制了计算的乘法深度。

    It is possible to choose the initial scale freely. One good strategy can be
    to is to set the initial scale S and primes P_i in the coeff_modulus to be
    very close to each other. If ciphertexts have scale S before multiplication,
    they have scale S^2 after multiplication, and S^2/P_i after rescaling. If all
    P_i are close to S, then S^2/P_i is close to S again. This way we stabilize the
    scales to be close to S throughout the computation. Generally, for a circuit
    of depth D, we need to rescale D times, i.e., we need to be able to remove D
    primes from the coefficient modulus. Once we have only one prime left in the
    coeff_modulus, the remaining prime must be larger than S by a few bits to
    preserve the pre-decimal-point value of the plaintext.
    自由选择初始scale。一个好策略是将原始scale S与coeff_modulus中的素数P_i设置得十分靠近。
    如果在乘法之前，密文缩放了S倍，则密文在乘法之后将缩放S^2倍，再重缩放之后变成了S^2/P_i。
    如果所有的P_i都接近S，那么S^2/P_也将接近S。通过这种方式，我们在整个计算过程中稳定了规模
    尺度，使其接近于S。通常对于一个深度为D的电路，我们需要重缩放D次，即我们需要在系数模中移
    除D个素数。一旦coeff_modulus中只有一个素数的时候，该素数一定要比S大几位，以保留明文的小
    数点前的值。


    Therefore, a generally good strategy is to choose parameters for the CKKS
    scheme as follows:

        (1) Choose a 60-bit prime as the first prime in coeff_modulus. This will
            give the highest precision when decrypting;
            选择一个60位的素数作为coeff_modulus中的第一个素数，这样将会在解密时拥有最高的准确率。

        (2) Choose another 60-bit prime as the last element of coeff_modulus, as
            this will be used as the special prime and should be as large as the
            largest of the other primes;
            选择另外一个60位的素数作为coeff_modulus中的最后一个元素，这样将会用作特殊素数，应该与其他最大的素数一样大。

        (3) Choose the intermediate primes to be close to each other.
            选择中间素数，大小尽可能互相靠近。

    We use CoeffModulus::Create to generate primes of the appropriate size. Note
    that our coeff_modulus is 200 bits total, which is below the bound for our
    poly_modulus_degree: CoeffModulus::MaxBitCount(8192) returns 218.
    */
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));

    /*
    We choose the initial scale to be 2^40. At the last level, this leaves us
    60-40=20 bits of precision before the decimal point, and enough (roughly
    10-20 bits) of precision after the decimal point. Since our intermediate
    primes are 40 bits (in fact, they are very close to 2^40), we can achieve
    scale stabilization as described above.
    我们选择一个初始scale为2^40，在最后一个level中，将留给我小数点前20位的精确位数，和足够的（粗略为10-20位）
    小数点后的精确位数。由于我们中间素数为40位（十分接近2^40），我们可以实现一个上述描述的scale 稳定。
    */
    double scale = pow(2.0, 40);

    //参数，打印
    auto context = SEALContext::Create(parms);
    print_parameters(context);
    cout << endl;

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

    //定义时间变量
    chrono::microseconds time_encode_input(0);
    chrono::microseconds time_encode_filter(0);
    chrono::microseconds time_decode_result(0);
    chrono::microseconds time_encrypt_input(0);
    chrono::microseconds time_decrypt_result(0);
    chrono::microseconds time_add(0);
    chrono::microseconds time_add_sum(0);
    chrono::microseconds time_multiply(0);
    chrono::microseconds time_packing_input(0);
    chrono::microseconds time_packing_filter(0);
    chrono::microseconds time_read_input(0);
    chrono::microseconds time_transform_input(0);
    chrono::microseconds time_rescale(0);
    chrono::microseconds time_relinearize(0);
    chrono::microseconds time_rotation(0);
    chrono::microseconds time_rotation_sum(0);
    chrono::microseconds time_sum(0);

    //输入部分，到时候如果改输入的话，可以修改这个部分

    ifstream f;
    f.open("C:\\Users\\DELL\\Desktop\\data.txt");
    string str;

    // 计算从文件中读取数据需要的时间
    time_start = chrono::high_resolution_clock::now();
    //二维输入数组
    vector<vector<double> > num;
    while (getline(f, str))
    {
        istringstream input(str);
        vector<double> tmp;
        double a;
        while (input >> a)
            tmp.push_back(a);

        num.push_back(tmp);
    }
    time_end = chrono::high_resolution_clock::now();
    time_read_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Read data from file:";
    cout << time_read_input.count() << " microseconds" << endl;

    // cout << endl;
    // print_line(__LINE__);
    // cout << endl;

    //转换后的输入矩阵
    vector<vector<double> > input_matrix;
    int count = 0;
    // 计算转换输入矩阵需要的时间
    time_start = chrono::high_resolution_clock::now();
    for (size_t h = 0; h < num.size() - 2 + 1; ++h)
    {
        for (size_t w = 0; w < num[0].size() - 2 + 1; ++w)
        {
            vector<double> tmp;
            for (size_t i = h; i < h + 2; ++i)
            {
                for (size_t j = w; j < w + 2; ++j)
                {
                    // cout << num[i][j] << " ";
                    tmp.push_back(num[i][j]);
                }
            }
            input_matrix.push_back(tmp);
            // cout << endl;
            count++;
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_transform_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Transform the input matrix: Done [";
    cout << time_transform_input.count() << " microseconds]" << endl;
    // cout << "共有" << count << "行" << endl;

    // 转换后的filter
    // 将filter packing到一条明文内
    /* cout << endl;
     print_line(__LINE__);
     cout << "转换后的filter vector:" << endl;*/
    double filter[4] = { 0.0, 1.0, 1.0, 0.0 };
    vector<double> filter_vector;

    // 计算将filter转换成对应的明文需要的时间
    time_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < count; ++j)
        {
            // cout << filter[i] << " ";

            filter_vector.push_back(filter[i]);
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_transform_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Transform the filter: Done [";
    cout << time_transform_input.count() << " microseconds]" << endl;

    //现将输入矩阵packing到一条密文内
    vector<double> input;
    input.reserve(slot_count);

    /*cout << endl << endl;
    print_line(__LINE__);
    cout << "处理后用于计算的input matrix:" << endl;
    cout << "[ ";*/

    // 计算将input packing起来需要的时间
    time_start = chrono::high_resolution_clock::now();
    //这里需要注意，是先将第一列列放进去，再放第二列，这样才能保证旋转后加和的点积的形式
    for (size_t i = 0; i < input_matrix[0].size(); ++i)
    {
        for (size_t j = 0; j < input_matrix.size(); ++j)
        {
            // cout << input_matrix[j][i] << ", ";
            input.push_back(input_matrix[j][i]);
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_packing_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    /*cout << " ]" << endl << endl;*/
    cout << "Packing input: Done [";
    cout << time_packing_input.count() << " microseconds]" << endl;

    // print_line(__LINE__);
    /*cout << "Input vector: " << endl;
    print_vector(input, 5, 7);*/

    cout << endl << "Evaluating input matrix and filter .........." << endl;

    /*
    We create plaintexts for PI, 0.4, and 1 using an overload of CKKSEncoder::encode
    that encodes the given floating-point value to every slot in the vector.
    我们使用CKKSEncoder::encode重载创建PI、0.4和1的明文，该重载将给定的浮点值编码到向量中的每个槽。
    */

    // 将filter vector进行编码，并计算编码的时间
    Plaintext filter_plain;
    time_start = chrono::high_resolution_clock::now();
    encoder.encode(filter_vector, scale, filter_plain);
    time_end = chrono::high_resolution_clock::now();
    time_encode_filter = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    // print_line(__LINE__);
    cout << "Encode filter: Done [";
    cout << time_encode_filter.count() << " microseconds]" << endl;

    Plaintext x_plain;
    // print_line(__LINE__);
    cout << "Encode input vectors：Done [";
    time_start = chrono::high_resolution_clock::now();
    encoder.encode(input, scale, x_plain); //将输入进行编码
    time_end = chrono::high_resolution_clock::now();
    time_encode_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_encode_input.count() << " microseconds]" << endl;

    Ciphertext x1_encrypted;
    // print_line(__LINE__);
    cout << "Encrypt input vectors：Done [";
    time_start = chrono::high_resolution_clock::now();
    encryptor.encrypt(x_plain, x1_encrypted); //将编码后的输入进行加密
    time_end = chrono::high_resolution_clock::now();
    time_encrypt_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_encrypt_input.count() << " microseconds]" << endl;

    /*
    To compute x^3 we first compute x^2 and relinearize. However, the scale has
    now grown to 2^80.重复线性化？
    */

    /*
    Now rescale; in addition to a modulus switch, the scale is reduced down by
    a factor equal to the prime that was switched away (40-bit prime). Hence, the
    new scale should be close to 2^40. Note, however, that the scale is not equal
    to 2^40: this is because the 40-bit prime is only close to 2^40.
    重缩放；除了模数转换，scale还会减少一个因子，这个因子等于被转换掉的素数(40位素数)。
    因此，新的scale应接近2^40。但是，请注意，这个scale并不等于2^40:这是因为这个40位的质数只接近2^40。
    */

    /*
    Now x3_encrypted is at a different level than x1_encrypted, which prevents us
    from multiplying them to compute x^3. We could simply switch x1_encrypted to
    the next parameters in the modulus switching chain. However, since we still
    need to multiply the x^3 term with PI (plain_coeff3), we instead compute PI*x
    first and multiply that with x^2 to obtain PI*x^3. To this end, we compute
    PI*x and rescale it back from scale 2^80 to something close to 2^40.
    这里还只计算到x^2，对于3次方，由于还需要乘以PI，因此先将x乘以PI。
    然后将结果重缩放到40位之后，再与x^2进行相乘。之后再进行relinearize以及重缩放。
    */
    // print_line(__LINE__);
    // cout << "Compute and rescale input*filter." << endl;
    // Ciphertext x1_encrypted_filter;
    // evaluator.multiply_plain(x1_encrypted, filter_plain, x1_encrypted_filter); // 加密的输入与未加密的filter进行相乘
    // cout << "    + Scale of input*filter before rescale: " << log2(x1_encrypted_filter.scale()) << " bits" << endl;
    // evaluator.rescale_to_next_inplace(x1_encrypted_filter);
    // cout << "    + Scale of input*filter after rescale: " << log2(x1_encrypted_filter.scale()) << " bits" << endl;
    // print_line(__LINE__);
    // cout << "Compute and rescale input*filter." << endl;
    Ciphertext x1_encrypted_filter;
    Plaintext plain;
    vector<double> result;

    // print_line(__LINE__);
    cout << "Input ciphertext multiplies filter plaintext：Done [";
    time_start = chrono::high_resolution_clock::now();
    evaluator.multiply_plain(x1_encrypted, filter_plain, x1_encrypted_filter); // 加密的输入与未加密的filter进行相乘
    time_end = chrono::high_resolution_clock::now();
    time_multiply = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_multiply.count() << " microseconds]" << endl;

    // print_line(__LINE__);
    cout << "Relinearize the result：Done [";
    time_start = chrono::high_resolution_clock::now();
    evaluator.relinearize_inplace(x1_encrypted_filter, relin_keys); // Relinearize 𝑟𝑒𝑠𝑢𝑙𝑡, keep same size
    time_end = chrono::high_resolution_clock::now();
    time_relinearize = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_relinearize.count() << " microseconds]" << endl << endl;

    /*print_line(__LINE__);
    cout << "加密的输入和未加密的filter相乘的结果：" << endl;
    decryptor.decrypt(x1_encrypted_filter, plain);

    encoder.decode(plain, result);
    cout << "[ ";
    for (int i = 0; i < result.size(); ++i)
    {
        printf("%.7f", result[i]);
        cout << ", ";
    }
    cout << " ]" << endl<<endl;

    encryptor.encrypt(plain, x1_encrypted_filter);
    */

    // Rotation
    Ciphertext rotated;

    // rotation，并将每次的结果与最终的结果输出。
    int c = 0;
    // print_line(__LINE__);
    // cout << "Rotation：";
    time_start = chrono::high_resolution_clock::now();
    for (int i = filter_vector.size() / 2; i >= count;)
    {
        c++;
        // print_line(__LINE__);
        if (i % 2 != 0)
        {
            i = int(i) + 1;
        }
        cout << "Half Rotation left " << c << ": " << i << " steps: Done [";
        time_start = chrono::high_resolution_clock::now();
        evaluator.rotate_vector(x1_encrypted_filter, i, gal_keys, rotated);
        time_end = chrono::high_resolution_clock::now();
        time_rotation = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_rotation.count() << " microseconds]" << endl;
        time_rotation_sum += time_rotation;

        cout << "Add: Done [";
        time_start = chrono::high_resolution_clock::now();
        evaluator.add_inplace(x1_encrypted_filter, rotated); //两者相加，结果存到x1_encrypted_filter中
        time_end = chrono::high_resolution_clock::now();
        time_add = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << time_add.count() << " microseconds]" << endl;
        time_add_sum += time_add;
        /*cout << "    + Decrypt and decode ...... Correct." << endl;*/

        /*decryptor.decrypt(x1_encrypted_filter, plain);
        encoder.decode(plain, result);
        cout << "[ ";
        for (int j = 0; j < i; ++j)
        {
            printf("%.7f", result[j]);
            cout << ", ";
        }
        cout << " ]" << endl << endl;
        encryptor.encrypt(plain, x1_encrypted_filter);*/
        i = i / 2;
    }
    // print_line(__LINE__);
    cout << endl << "Sum of rotation time: [";
    cout << time_rotation_sum.count() << " microseconds]" << endl;
    // print_line(__LINE__);
    cout << "Sum of add time: [";
    cout << time_add_sum.count() << " microseconds]" << endl << endl;

    /*
    Since x3_encrypted and x1_encrypted_coeff3 have the same exact scale and use
    the same encryption parameters, we can multiply them together. We write the
    result to x3_encrypted, relinearize, and rescale. Note that again the scale
    is something close to 2^40, but not exactly 2^40 due to yet another scaling
    by a prime. We are down to the last level in the modulus switching chain.
    */

    /*
    Next we compute the degree one term. All this requires is one multiply_plain
    with plain_coeff1. We overwrite x1_encrypted with the result.
    计算一次方项。所有这些都需要一个multiply_plain和plain_coeff1。我们用结果覆盖x1_encrypted。
    */

    /*
    Now we would hope to compute the sum of all three terms. However, there is
    a serious problem: the encryption parameters used by all three terms are
    different due to modulus switching from rescaling.
    将三项加起来，但是有一个问题：由于从缩放到模转换，所以这三项所使用的加密参数是不同的。

    Encrypted addition and subtraction require that the scales of the inputs are
    the same, and also that the encryption parameters (parms_id) match. If there
    is a mismatch, Evaluator will throw an exception.
    加密加法和减法要求输入的scale是相同的，同时加密参数也要匹配。
    如果不匹配，Evaluator将会抛出一个异常。
    */
    /* cout << endl;
     print_line(__LINE__);
     cout << "Parameters used by all three terms are different." << endl;
     cout << "    + Modulus chain index for x3_encrypted: "
          << context->get_context_data(x3_encrypted.parms_id())->chain_index() << endl;
     cout << "    + Modulus chain index for x1_encrypted: "
          << context->get_context_data(x1_encrypted.parms_id())->chain_index() << endl;
     cout << "    + Modulus chain index for plain_coeff0: "
          << context->get_context_data(plain_coeff0.parms_id())->chain_index() << endl;
     cout << endl;*/

    /*
    Let us carefully consider what the scales are at this point. We denote the
    primes in coeff_modulus as P_0, P_1, P_2, P_3, in this order. P_3 is used as
    the special modulus and is not involved in rescalings. After the computations
    above the scales in ciphertexts are:

        - Product x^2 has scale 2^80 and is at level 2;
        - Product PI*x has scale 2^80 and is at level 2;
        - We rescaled both down to scale 2^80/P_2 and level 1;
        - Product PI*x^3 has scale (2^80/P_2)^2;
        - We rescaled it down to scale (2^80/P_2)^2/P_1 and level 0;
        - Product 0.4*x has scale 2^80;
        - We rescaled it down to scale 2^80/P_2 and level 1;
        - The contant term 1 has scale 2^40 and is at level 2.

    Although the scales of all three terms are approximately 2^40, their exact
    values are different, hence they cannot be added together.
    */
    /* print_line(__LINE__);
     cout << "The exact scales of all three terms are different:" << endl;
     ios old_fmt(nullptr);
     old_fmt.copyfmt(cout);
     cout << fixed << setprecision(10);
     cout << "    + Exact scale in PI*x^3: " << x3_encrypted.scale() << endl;
     cout << "    + Exact scale in  0.4*x: " << x1_encrypted.scale() << endl;
     cout << "    + Exact scale in      1: " << plain_coeff0.scale() << endl;
     cout << endl;
     cout.copyfmt(old_fmt);*/

    /*
    There are many ways to fix this problem. Since P_2 and P_1 are really close
    to 2^40, we can simply "lie" to Microsoft SEAL and set the scales to be the
    same. For example, changing the scale of PI*x^3 to 2^40 simply means that we
    scale the value of PI*x^3 by 2^120/(P_2^2*P_1), which is very close to 1.
    This should not result in any noticeable error.

    Another option would be to encode 1 with scale 2^80/P_2, do a multiply_plain
    with 0.4*x, and finally rescale. In this case we would need to additionally
    make sure to encode 1 with appropriate encryption parameters (parms_id).

    In this example we will use the first (simplest) approach and simply change
    the scale of PI*x^3 and 0.4*x to 2^40.
    */
    /*print_line(__LINE__);
    cout << "Normalize scales to 2^40." << endl;
    x3_encrypted.scale() = pow(2.0, 40);
    x1_encrypted.scale() = pow(2.0, 40);*/

    /*
    We still have a problem with mismatching encryption parameters. This is easy
    to fix by using traditional modulus switching (no rescaling). CKKS supports
    modulus switching just like the BFV scheme, allowing us to switch away parts
    of the coefficient modulus when it is simply not needed.
    CKKS支持模量切换，就像BFV方案一样，允许我们在不需要的时候切换部分系数模量。
    */
    // print_line(__LINE__);
    // cout << "Normalize encryption parameters to the lowest level." << endl;
    ////将所有的都转换为一样的模数系数id了
    // parms_id_type last_parms_id = x3_encrypted.parms_id();
    // evaluator.mod_switch_to_inplace(x1_encrypted, last_parms_id);
    // evaluator.mod_switch_to_inplace(plain_coeff0, last_parms_id);

    ///*
    // All three ciphertexts are now compatible and can be added.
    //*/
    // print_line(__LINE__);
    // cout << "Compute PI*x^3 + 0.4*x + 1." << endl;
    // Ciphertext encrypted_result;
    // evaluator.add(x3_encrypted, x1_encrypted, encrypted_result);
    // evaluator.add_plain_inplace(encrypted_result, plain_coeff0);

    /*
    First print the true result.
    */

    /*print_line(__LINE__);
    cout << "Decrypt and decode input*filter." << endl;

    */
    // print_line(__LINE__);
    cout << "Expected result:" << endl;
    vector<double> true_result;
    cout << "[ ";

    time_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < input_matrix.size(); i++)
    {
        double x = 0;
        for (size_t j = 0; j < input_matrix[0].size(); ++j)
        {
            x += input_matrix[i][j] * filter[j];
        }
        printf("%.7f", x);
        cout << ", ";

        true_result.push_back(x);
    }
    cout << " ]" << endl;
    // print_line(__LINE__);
    time_end = chrono::high_resolution_clock::now();
    time_transform_input = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    cout << "Generating expected result: Dome [" << time_transform_input.count() << " microseconds]" << endl;

    // print_line(__LINE__);

    /*
    Decrypt, decode, and print the result.
    */

    Plaintext plain_result;
    // print_line(__LINE__);
    cout << endl << "Decrypt the result: Done [";
    time_start = chrono::high_resolution_clock::now();
    decryptor.decrypt(x1_encrypted_filter, plain_result);
    time_end = chrono::high_resolution_clock::now();
    time_decrypt_result = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_decrypt_result.count() << " microseconds]" << endl;

    // print_line(__LINE__);
    cout << "Decode the result: Done [";
    time_start = chrono::high_resolution_clock::now();
    encoder.decode(plain_result, result);
    time_end = chrono::high_resolution_clock::now();
    time_decode_result = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << time_decode_result.count() << " microseconds]" << endl;

    cout << endl << "Computed result ...... Correct." << endl;
    cout << "[ ";
    for (int i = 0; i < input_matrix.size(); ++i)
    {
        printf("%.7f", result[i]);
        cout << ", ";
    }
    cout << " ]" << endl << endl;

    cout << "Sum of encrypted input * filter time: [";
    time_sum = time_packing_input + time_encode_filter + time_encode_input + time_add_sum + time_decode_result +
               time_decrypt_result + time_encrypt_input + time_multiply + time_encrypt_input + time_packing_filter +
               time_relinearize + time_rotation_sum;
    cout << time_sum.count() << " microseconds]" << endl << endl;
    /*
    While we did not show any computations on complex numbers in these examples,
    the CKKSEncoder would allow us to have done that just as easily. Additions
    and multiplications of complex numbers behave just as one would expect.
    */
}
