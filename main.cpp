#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// 定义像片坐标类
class Photo {
public:
    double u, v;
    Photo(double u, double v) : u(u), v(v) {}
};

// 像片坐标的单位转换函数
Photo create_photo(double x, double y) {
    return Photo(x / 1000, y / 1000);
}

// 定义地面点类
class Landspace {
public:
    double x, y, z;
    Landspace(double x, double y, double z) : x(x), y(y), z(z) {}

    // 计算地面点坐标x的平均值
    static double avg_x(const vector<Landspace>& lands) {
        double sum = 0;
        for (const auto& land : lands) sum += land.x;
        return sum / lands.size();
    }

    // 计算地面点坐标y的平均值
    static double avg_y(const vector<Landspace>& lands) {
        double sum = 0;
        for (const auto& land : lands) sum += land.y;
        return sum / lands.size();
    }

    // 计算地面点坐标z的平均值
    static double avg_z(const vector<Landspace>& lands) {
        double sum = 0;
        for (const auto& land : lands) sum += land.z;
        return sum / lands.size();
    }
};

// 定义旋转矩阵类
class RMatrix {
public:
    Matrix3d R;
    // 构造函数，根据phi、omega、kappa计算旋转矩阵
    RMatrix(double phi, double omega, double kappa) {
        // 计算phi、omega、kappa的cos和sin值
        double cos_phi = cos(phi), sin_phi = sin(phi);
        double cos_omega = cos(omega), sin_omega = sin(omega);
        double cos_kappa = cos(kappa), sin_kappa = sin(kappa);

        // 计算旋转矩阵
        R(0, 0) = cos_phi * cos_kappa - sin_phi * sin_omega * sin_kappa;
        R(0, 1) = -cos_phi * sin_kappa - sin_phi * sin_omega * cos_kappa;
        R(0, 2) = -sin_phi * cos_omega;
        R(1, 0) = cos_omega * sin_kappa;
        R(1, 1) = cos_omega * cos_kappa;
        R(1, 2) = -sin_omega;
        R(2, 0) = sin_phi * cos_kappa + cos_phi * sin_omega * sin_kappa;
        R(2, 1) = -sin_phi * sin_kappa + cos_phi * sin_omega * cos_kappa;
        R(2, 2) = cos_phi * cos_omega;
    }
};

// 定义A矩阵类
class AMatrix {
public:
    Matrix<double, 2, 6> A;
    // 根据像片坐标和R矩阵、焦距等参数计算A矩阵
    void create_AMatrix(double x, double y, const RMatrix& R, double _X_, double _Y_, double _Z_, double f, double x0, double y0, double phi, double omega, double kappa) {
        double dx = x - x0, dy = y - y0;
        A(0, 0) = (R.R(0, 0) * f + R.R(0, 2) * dx) / _Z_;
        A(0, 1) = (R.R(1, 0) * f + R.R(1, 2) * dx) / _Z_;
        A(0, 2) = (R.R(2, 0) * f + R.R(2, 2) * dx) / _Z_;
        A(1, 0) = (R.R(0, 1) * f + R.R(0, 2) * dy) / _Z_;
        A(1, 1) = (R.R(1, 1) * f + R.R(1, 2) * dy) / _Z_;
        A(1, 2) = (R.R(2, 1) * f + R.R(2, 2) * dy) / _Z_;

        // 简写三角函数
        double cos_kappa = cos(kappa), sin_kappa = sin(kappa);
        double cos_omega = cos(omega), sin_omega = sin(omega);

        // 计算A矩阵的后三列
        A(0, 3) = dy * sin_omega - (((dx * cos_kappa - dy * sin_kappa) * dx / f) + f * cos_kappa) * cos_omega;
        A(0, 4) = -f * sin_kappa - dx / f * (dx * sin_kappa + dy * cos_kappa);
        A(0, 5) = dy;
        A(1, 3) = -dx * sin_omega - (((dx * cos_kappa - dy * sin_kappa) * dy / f) - f * sin_kappa) * cos_omega;
        A(1, 4) = -f * cos_kappa - dy / f * (dx * sin_kappa + dy * cos_kappa);
        A(1, 5) = -dx;
    }
};

/**
 * @brief 计算旋转矩阵
 *
 * 该函数根据给定的外方位角phi、omega和kappa计算旋转矩阵R。
 *
 * @param phi 外方位角phi（弧度）
 * @param omega 外方位角omega（弧度）
 * @param kappa 外方位角kappa（弧度）
 * @param R 计算得到的旋转矩阵
 *
 * 该函数的主要步骤如下：
 * 1. 使用给定的外方位角phi、omega和kappa初始化旋转矩阵R。
 * 2. 旋转矩阵R的计算基于三角函数的组合。
 *
 * @note 该函数假设输入的外方位角phi、omega和kappa均为弧度制。
 *
 * @warning 如果输入的外方位角不是弧度制，计算结果将不正确。
 */
void calculate_rotation_matrix(double phi, double omega, double kappa, RMatrix& R) {
    cout << "calculate_rotation_matrix..." << endl;

    R = RMatrix(phi, omega, kappa);

    cout << "Rotation matrix: " << endl;
    cout << R.R << endl;
}

/**
 * @brief 计算近似值
 *
 * 该函数计算每个照片点的近似值，包括旋转矩阵R的应用和投影坐标的计算。
 *
 * @param photos 包含照片点的向量，每个点包含u和v坐标
 * @param lands 包含地面点的向量，每个点包含x, y, z坐标
 * @param Xs 摄影中心的X坐标
 * @param Ys 摄影中心的Y坐标
 * @param Zs 摄影中心的Z坐标
 * @param f 摄影机焦距
 * @param x0 摄影机主点的X坐标
 * @param y0 摄影机主点的Y坐标
 * @param R 旋转矩阵
 * @param Xclose 计算得到的近似X坐标的向量
 * @param Yclose 计算得到的近似Y坐标的向量
 * @param X 计算得到的X坐标的向量
 * @param Y 计算得到的Y坐标的向量
 * @param Z 计算得到的Z坐标的向量
 *
 * 该函数的主要步骤如下：
 * 1. 对于每个照片点，计算其在摄影中心坐标系下的坐标X, Y, Z。
 * 2. 使用旋转矩阵R将地面点坐标转换到摄影中心坐标系。
 * 3. 计算投影坐标Xclose和Yclose。
 * 4. 如果Z坐标不为零，则计算投影坐标，否则将投影坐标设置为主点坐标。
 *
 * @note 该函数假设输入向量photos和lands的大小相同。
 *
 * @warning 如果Z坐标为零，投影坐标将被设置为主点坐标，这可能导致不准确的结果。
 */
void calculate_approximate_values(const vector<Photo>& photos, const vector<Landspace>& lands, double Xs, double Ys, double Zs, double f, double x0, double y0, const RMatrix& R, vector<double>& Xclose, vector<double>& Yclose, vector<double>& X, vector<double>& Y, vector<double>& Z) {
    cout << "calculate_approximate_values..." << endl;
    for (size_t i = 0; i < photos.size(); ++i) {
        X[i] = R.R(0, 0) * (lands[i].x - Xs) + R.R(1, 0) * (lands[i].y - Ys) + R.R(2, 0) * (lands[i].z - Zs);
        Y[i] = R.R(0, 1) * (lands[i].x - Xs) + R.R(1, 1) * (lands[i].y - Ys) + R.R(2, 1) * (lands[i].z - Zs);
        Z[i] = R.R(0, 2) * (lands[i].x - Xs) + R.R(1, 2) * (lands[i].y - Ys) + R.R(2, 2) * (lands[i].z - Zs);

        if (Z[i] != 0) {
            Xclose[i] = x0 - f * X[i] / Z[i];
            Yclose[i] = y0 - f * Y[i] / Z[i];
        } else {
            Xclose[i] = x0;
            Yclose[i] = y0;
        }

        cout << "Photo " << i << ": " << X[i] << ", " << Y[i] << ", " << Z[i] << endl;
        cout << "Photo " << i << " close: " << Xclose[i] << ", " << Yclose[i] << endl;
    }
}

/**
 * @brief 计算误差方程系数
 *
 * 该函数计算每个照片点的误差方程系数，包括旋转矩阵R的应用和误差向量L的计算。
 *
 * @param photos 包含照片点的向量，每个点包含u和v坐标
 * @param lands 包含地面点的向量，每个点包含x, y, z坐标
 * @param Xclose 计算得到的近似X坐标的向量
 * @param Yclose 计算得到的近似Y坐标的向量
 * @param X 计算得到的X坐标的向量
 * @param Y 计算得到的Y坐标的向量
 * @param Z 计算得到的Z坐标的向量
 * @param f 摄影机焦距
 * @param x0 摄影机主点的X坐标
 * @param y0 摄影机主点的Y坐标
 * @param phi 外方位角phi
 * @param omega 外方位角omega
 * @param kappa 外方位角kappa
 * @param R 旋转矩阵
 * @param A_matrices 计算得到的误差方程系数矩阵的向量
 * @param L 计算得到的误差向量
 *
 * 该函数的主要步骤如下：
 * 1. 对于每个照片点，计算其误差方程系数矩阵A。
 * 2. 使用旋转矩阵R和地面点坐标计算误差方程系数。
 * 3. 计算误差向量L。
 *
 * @note 该函数假设输入向量photos和lands的大小相同。
 *
 * @warning 如果输入数据不一致，可能会导致计算错误。
 */
void calculate_error_equation_coefficients(const vector<Photo>& photos, const vector<Landspace>& lands, const vector<double>& Xclose, const vector<double>& Yclose, const vector<double>& X, const vector<double>& Y, const vector<double>& Z, double f, double x0, double y0, double phi, double omega, double kappa, const RMatrix& R, vector<AMatrix>& A_matrices, VectorXd& L) {
    cout << "calculate_error_equation_coefficients..." << endl;

    for (size_t i = 0; i < photos.size(); ++i) {
        A_matrices[i].create_AMatrix(Xclose[i], Yclose[i], R, X[i], Y[i], Z[i], f, x0, y0, phi, omega, kappa);
        L(2 * i) = photos[i].u - Xclose[i];
        L(2 * i + 1) = photos[i].v - Yclose[i];

        cout << "Photo " << i << ": " << L(2 * i) << ", " << L(2 * i + 1) << endl;
    }
}

/**
 * @brief 计算外方位元素修正值
 *
 * 该函数计算外方位元素的修正值，使用最小二乘法求解误差方程。
 *
 * @param A_matrices 包含误差方程系数矩阵的向量
 * @param L 误差向量
 * @param X 计算得到的外方位元素修正值向量
 *
 * 该函数的主要步骤如下：
 * 1. 构建误差方程系数矩阵A。
 * 2. 使用最小二乘法求解AX = L，得到外方位元素修正值X。
 *
 * @note 该函数假设输入向量A_matrices和L的大小一致。
 *
 * @warning 如果输入数据不一致，可能会导致计算错误。
 */
void calculate_exterior_orientation_corrections(const vector<AMatrix>& A_matrices, const VectorXd& L, VectorXd& X) {
    cout << "calculate_exterior_orientation_corrections..." << endl;

    size_t photo_count = A_matrices.size();
    MatrixXd A(2 * photo_count, 6);

    for (size_t i = 0; i < photo_count; ++i) {
        A.block<2, 6>(2 * i, 0) = A_matrices[i].A;
    }

    X = (A.transpose() * A).ldlt().solve(A.transpose() * L);

    cout << "Exterior orientation corrections: " << X << endl;
}

/**
 * @brief 检查是否收敛
 *
 * 该函数检查外方位元素修正值是否收敛，判断条件是修正值的绝对值的最大值是否小于给定的容差。
 *
 * @param X 外方位元素修正值向量
 * @param tolerance 收敛容差
 * @return 如果所有修正值的绝对值都小于容差，则返回true；否则返回false
 *
 * 该函数的主要步骤如下：
 * 1. 计算修正值向量X中每个元素的绝对值。
 * 2. 找到绝对值的最大值。
 * 3. 判断最大值是否小于给定的容差。
 *
 * @note 该函数假设输入向量X的大小为6。
 *
 * @warning 如果输入向量X的大小不为6，可能会导致计算错误。
 */
bool check_convergence(const VectorXd& X, double tolerance) {
    cout << "check_convergence..." << endl;

    return X.cwiseAbs().maxCoeff() < tolerance;

    cout << "Convergence: " << (X.cwiseAbs().maxCoeff() < tolerance) << endl;
}

/**
 * @brief 计算单位权中误差
 *
 * 该函数计算单位权中误差（unit weight standard deviation），用于评估外方位元素修正值的精度。
 *
 * @param A_matrices 包含误差方程系数矩阵的向量
 * @param L 误差向量
 * @param X 计算得到的外方位元素修正值向量
 * @return 单位权中误差
 *
 * 该函数的主要步骤如下：
 * 1. 初始化残差向量V为零向量。
 * 2. 对于每个照片点，计算残差向量V。
 * 3. 计算残差向量V的平方和。
 * 4. 计算单位权中误差。
 *
 * @note 该函数假设输入向量A_matrices和L的大小一致。
 *
 * @warning 如果输入数据不一致，可能会导致计算错误。
 */
double calculate_unit_weight_standard_deviation(const vector<AMatrix>& A_matrices, const VectorXd& L, const VectorXd& X) {
    cout << "calculate_unit_weight_standard_deviation..." << endl;

    size_t photo_count = A_matrices.size();
    VectorXd V = VectorXd::Zero(2 * photo_count);

    for (size_t i = 0; i < photo_count; ++i) {
        V.segment<2>(2 * i) = A_matrices[i].A * X - L.segment<2>(2 * i);
    }

    double vv = V.squaredNorm();

    cout << "Residual vector: " << V.transpose() << endl;
    cout << "Sum of squares of residuals: " << vv << endl;
    return sqrt(vv / (2 * photo_count - 6));
}

int main() {
    // 初始化像片坐标和地面点坐标
    vector<Photo> photos = {
        create_photo(-86.15, -68.99),
        create_photo(-53.40, 82.21),
        create_photo(-14.78, -76.63),
        create_photo(10.46, 64.43)
    };

    vector<Landspace> lands = {
        Landspace(36589.41, 25273.32, 2195.17),
        Landspace(37631.08, 31324.51, 728.69),
        Landspace(39100.97, 24934.98, 2386.50),
        Landspace(40426.54, 30319.81, 757.31)
    };

    // 初始化参数
    const double f = 153.24 / 1000;
    const double x0 = 0;
    const double y0 = 0;
    const double m = 50000;

    double Xs = Landspace::avg_x(lands);
    double Ys = Landspace::avg_y(lands);
    double Zs = f * m + Landspace::avg_z(lands);
    double phi0 = 0;
    double omega0 = 0;
    double kappa0 = 0;

    // 初始化旋转矩阵
    RMatrix R(phi0, omega0, kappa0);

    // 初始化近似值
    vector<double> Xclose(photos.size());
    vector<double> Yclose(photos.size());
    vector<double> _X_(photos.size());
    vector<double> _Y_(photos.size());
    vector<double> _Z_(photos.size());

    // 初始化误差方程系数
    calculate_approximate_values(photos, lands, Xs, Ys, Zs, f, x0, y0, R, Xclose, Yclose, _X_, _Y_, _Z_);

    // 初始化A矩阵和L矩阵
    vector<AMatrix> A_matrices(photos.size());
    VectorXd L(2 * photos.size());

    // 初始化外方位元素修正值
    calculate_error_equation_coefficients(photos, lands, Xclose, Yclose, _X_, _Y_, _Z_, f, x0, y0, phi0, omega0, kappa0, R, A_matrices, L);

    // 计算外方位元素修正值
    VectorXd X(6);

    // 迭代计算外方位元素修正值
    calculate_exterior_orientation_corrections(A_matrices, L, X);

    Xs += X[0];
    Ys += X[1];
    Zs += X[2];
    phi0 += X[3];
    omega0 += X[4];
    kappa0 += X[5];

    double tolerance = 1 * M_PI / 180 / 60;
    int max_iterations = 10;
    int iteration = 0;

    // 迭代计算外方位元素修正值
    while (iteration < max_iterations) {
        calculate_rotation_matrix(phi0, omega0, kappa0, R);
        calculate_approximate_values(photos, lands, Xs, Ys, Zs, f, x0, y0, R, Xclose, Yclose, _X_, _Y_, _Z_);
        calculate_error_equation_coefficients(photos, lands, Xclose, Yclose, _X_, _Y_, _Z_, f, x0, y0, phi0, omega0, kappa0, R, A_matrices, L);
        calculate_exterior_orientation_corrections(A_matrices, L, X);

        Xs += X[0];
        Ys += X[1];
        Zs += X[2];
        phi0 += X[3];
        omega0 += X[4];
        kappa0 += X[5];

        cout << "Iteration " << iteration << ": " << X.transpose() << endl;
        cout << "Xs: " << Xs << ", Ys: " << Ys << ", Zs: " << Zs << ", phi0: " << phi0 << ", omega0: " << omega0 << ", kappa0: " << kappa0 << endl;

        if (check_convergence(X, tolerance)) {
            break;
        }

        iteration++;
    }

    if (iteration >= max_iterations) {
        cout << "Failed to converge after " << max_iterations << " iterations." << endl;
    }

    // 计算单位权中误差
    double m0 = calculate_unit_weight_standard_deviation(A_matrices, L, X);
    cout << "Unit weight standard deviation (m0): " << m0 << endl;
    cout << R.R << endl;

    return 0;
}