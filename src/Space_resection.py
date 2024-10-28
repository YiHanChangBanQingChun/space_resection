import numpy as np

# 定义像片坐标
photos = [
    (-86.15 / 1000, -68.99 / 1000),
    (-53.40 / 1000, 82.21 / 1000),
    (-14.78 / 1000, -76.63 / 1000),
    (10.46 / 1000, 64.43 / 1000)
]

# 定义地面点坐标
lands = [
    (36589.41, 25273.32, 2195.17),
    (37631.08, 31324.51, 728.69),
    (39100.97, 24934.98, 2386.50),
    (40426.54, 30319.81, 757.31)
]

# 计算地面点坐标的平均值
def avg_coord(coords, index):
    return sum(coord[index] for coord in coords) / len(coords)

# 初始化参数
f = 153.24 / 1000
x0 = 0
y0 = 0
m = 50000

Xs = avg_coord(lands, 0)
Ys = avg_coord(lands, 1)
Zs = f * m + avg_coord(lands, 2)
phi0 = 0
omega0 = 0
kappa0 = 0

# 计算旋转矩阵
def calculate_rotation_matrix(phi, omega, kappa):
    cos_phi, sin_phi = np.cos(phi), np.sin(phi)
    cos_omega, sin_omega = np.cos(omega), np.sin(omega)
    cos_kappa, sin_kappa = np.cos(kappa), np.sin(kappa)

    R = np.array([
        [cos_phi * cos_kappa - sin_phi * sin_omega * sin_kappa, -cos_phi * sin_kappa - sin_phi * sin_omega * cos_kappa, -sin_phi * cos_omega],
        [cos_omega * sin_kappa, cos_omega * cos_kappa, -sin_omega],
        [sin_phi * cos_kappa + cos_phi * sin_omega * sin_kappa, -sin_phi * sin_kappa + cos_phi * sin_omega * cos_kappa, cos_phi * cos_omega]
    ])
    return R

# 计算近似值
def calculate_approximate_values(photos, lands, Xs, Ys, Zs, f, x0, y0, R):
    Xclose, Yclose, X, Y, Z = [], [], [], [], []
    for i in range(len(photos)):
        X_i = R[0, 0] * (lands[i][0] - Xs) + R[1, 0] * (lands[i][1] - Ys) + R[2, 0] * (lands[i][2] - Zs)
        Y_i = R[0, 1] * (lands[i][0] - Xs) + R[1, 1] * (lands[i][1] - Ys) + R[2, 1] * (lands[i][2] - Zs)
        Z_i = R[0, 2] * (lands[i][0] - Xs) + R[1, 2] * (lands[i][1] - Ys) + R[2, 2] * (lands[i][2] - Zs)
        X.append(X_i)
        Y.append(Y_i)
        Z.append(Z_i)

        if Z_i != 0:
            Xclose.append(x0 - f * X_i / Z_i)
            Yclose.append(y0 - f * Y_i / Z_i)
        else:
            Xclose.append(x0)
            Yclose.append(y0)
    return Xclose, Yclose, X, Y, Z

# 计算误差方程系数
def calculate_error_equation_coefficients(photos, Xclose, Yclose, X, Y, Z, f, x0, y0, phi, omega, kappa, R):
    A_matrices = []
    L = []
    for i in range(len(photos)):
        dx = Xclose[i] - x0
        dy = Yclose[i] - y0
        A = np.zeros((2, 6))
        A[0, 0] = (R[0, 0] * f + R[0, 2] * dx) / Z[i]
        A[0, 1] = (R[1, 0] * f + R[1, 2] * dx) / Z[i]
        A[0, 2] = (R[2, 0] * f + R[2, 2] * dx) / Z[i]
        A[1, 0] = (R[0, 1] * f + R[0, 2] * dy) / Z[i]
        A[1, 1] = (R[1, 1] * f + R[1, 2] * dy) / Z[i]
        A[1, 2] = (R[2, 1] * f + R[2, 2] * dy) / Z[i]

        cos_kappa, sin_kappa = np.cos(kappa), np.sin(kappa)
        cos_omega, sin_omega = np.cos(omega), np.sin(omega)

        A[0, 3] = dy * sin_omega - (((dx * cos_kappa - dy * sin_kappa) * dx / f) + f * cos_kappa) * cos_omega
        A[0, 4] = -f * sin_kappa - dx / f * (dx * sin_kappa + dy * cos_kappa)
        A[0, 5] = dy
        A[1, 3] = -dx * sin_omega - (((dx * cos_kappa - dy * sin_kappa) * dy / f) - f * sin_kappa) * cos_omega
        A[1, 4] = -f * cos_kappa - dy / f * (dx * sin_kappa + dy * cos_kappa)
        A[1, 5] = -dx

        A_matrices.append(A)
        L.append(photos[i][0] - Xclose[i])
        L.append(photos[i][1] - Yclose[i])
    return A_matrices, np.array(L)

# 计算外方位元素修正值
def calculate_exterior_orientation_corrections(A_matrices, L):
    A = np.vstack(A_matrices)
    X = np.linalg.lstsq(A, L, rcond=None)[0]
    return X

# 检查收敛性
def check_convergence(X, tolerance):
    return np.max(np.abs(X)) < tolerance

# 计算单位权中误差
def calculate_unit_weight_standard_deviation(A_matrices, L, X):
    V = np.hstack([A @ X - L[i:i+2] for i, A in enumerate(A_matrices)])
    vv = np.sum(V**2)
    return np.sqrt(vv / (2 * len(A_matrices) - 6))

# 初始化旋转矩阵
R = calculate_rotation_matrix(phi0, omega0, kappa0)

# 初始化近似值
Xclose, Yclose, _X_, _Y_, _Z_ = calculate_approximate_values(photos, lands, Xs, Ys, Zs, f, x0, y0, R)

# 初始化误差方程系数
A_matrices, L = calculate_error_equation_coefficients(photos, Xclose, Yclose, _X_, _Y_, _Z_, f, x0, y0, phi0, omega0, kappa0, R)

# 计算外方位元素修正值
X = calculate_exterior_orientation_corrections(A_matrices, L)

Xs += X[0]
Ys += X[1]
Zs += X[2]
phi0 += X[3]
omega0 += X[4]
kappa0 += X[5]

tolerance = 1 * np.pi / 180 / 60
max_iterations = 10
iteration = 0

# 迭代计算外方位元素修正值
while iteration < max_iterations:
    R = calculate_rotation_matrix(phi0, omega0, kappa0)
    Xclose, Yclose, _X_, _Y_, _Z_ = calculate_approximate_values(photos, lands, Xs, Ys, Zs, f, x0, y0, R)
    A_matrices, L = calculate_error_equation_coefficients(photos, Xclose, Yclose, _X_, _Y_, _Z_, f, x0, y0, phi0, omega0, kappa0, R)
    X = calculate_exterior_orientation_corrections(A_matrices, L)

    Xs += X[0]
    Ys += X[1]
    Zs += X[2]
    phi0 += X[3]
    omega0 += X[4]
    kappa0 += X[5]

    print(f"Iteration {iteration}: {X}")
    print(f"Xs: {Xs}, Ys: {Ys}, Zs: {Zs}, phi0: {phi0}, omega0: {omega0}, kappa0: {kappa0}")

    if check_convergence(X, tolerance):
        break

    iteration += 1

if iteration >= max_iterations:
    print(f"Failed to converge after {max_iterations} iterations.")

# 计算单位权中误差
m0 = calculate_unit_weight_standard_deviation(A_matrices, L, X)
print(f"Unit weight standard deviation (m0): {m0}")
print(R)