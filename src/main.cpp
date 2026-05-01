#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// 全局物理常数
const double GAMMA = 1.4;
const double R = 1.01e5 / 1.23 / 286.1;

// 几何参数
const double E = 10.0;
const double THETA = 5.352 * M_PI / 180.0;

// 计算域
const double H = 40.0;
const double L = 65.0;

// 初始条件
const double Ma1 = 2.0;
const double p1 = 1.01e5;
const double T1 = 286.1;
const double rho1 = 1.23;

// 数值参数
const double C = 0.5;
const double Cy = 0.6;
const int maxIteration = 10000;

// 网格参数
const double max_eta = 1.0;
const int num_eta = 41;
const double d_eta = max_eta / (num_eta - 1);

vector<double> linspace(double start, double end, int num) {
    vector<double> v(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) v[i] = start + i * step;
    return v;
}

double PrandtlMeyerFunc(double Ma, double f2) {
    double s = sqrt((GAMMA + 1) / (GAMMA - 1));
    double t = sqrt((GAMMA - 1) / (GAMMA + 1) * (Ma * Ma - 1));
    return s * atan(t) - atan(sqrt(Ma * Ma - 1)) - f2;
}

struct AnalyticSolution {
    double Ma, p, T, rho, u, v;
};
AnalyticSolution solution_analytic(double Ma1, double p1, double T1, double rho1, double phi, double theta) {
    double f1 = PrandtlMeyerFunc(Ma1, 0.0);
    double f2 = f1 + phi;
    double xm[3] = {1.0, 2.5, 5.0};
    double eps = 1e-5;
    int n = 0;
    while (fabs(xm[0] - xm[2]) >= eps) {
        double fm[3];
        for (int i = 0; i < 3; ++i) fm[i] = PrandtlMeyerFunc(xm[i], f2);
        if (fm[1] * fm[2] < 0) {
            xm[0] = xm[1];
            xm[1] = (xm[1] + xm[2]) / 2.0;
        } else if (fm[1] * fm[2] > 0) {
            xm[2] = xm[1];
            xm[1] = (xm[0] + xm[1]) / 2.0;
        } else {
            break;
        }
        if (++n > 1000) {
            cerr << "Warning: Max iteration for Ma2" << endl;
            break;
        }
    }
    double Ma2 = xm[1];
    double p2 = p1 * pow((1 + (GAMMA - 1) / 2 * Ma1 * Ma1) / (1 + (GAMMA - 1) / 2 * Ma2 * Ma2), GAMMA / (GAMMA - 1));
    double T2 = T1 * (1 + (GAMMA - 1) / 2 * Ma1 * Ma1) / (1 + (GAMMA - 1) / 2 * Ma2 * Ma2);
    double rho2 = p2 / (R * T2);
    double u2 = Ma2 * sqrt(GAMMA * R * T2 / (1 + tan(theta) * tan(theta)));
    double v2 = -u2 * tan(theta);
    return {Ma2, p2, T2, rho2, u2, v2};
}

void height(double x, double eta, double& h, double& P_eta_x, double& angle) {
    if (x <= E) {
        h = H;
        P_eta_x = 0.0;
        angle = 0.0;
    } else {
        h = H + (x - E) * tan(THETA);
        P_eta_x = (1 - eta) * tan(THETA) / h;
        angle = THETA;
    }
}

vector<double> calculate_F(double rho, double u, double v, double p) {
    vector<double> F(4);
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = GAMMA / (GAMMA - 1) * p * u + rho * u * (u * u + v * v) / 2.0;
    return F;
}

vector<double> calculate_G(const vector<double>& F, double rho) {
    double F1 = F[0], F2 = F[1], F3 = F[2];
    vector<double> G(4);
    G[0] = rho * F3 / F1;
    G[1] = F3;
    G[2] = rho * (F3 / F1) * (F3 / F1) + F2 - F1 * F1 / rho;
    G[3] = GAMMA / (GAMMA - 1) * (F2 - F1 * F1 / rho) * (F3 / F1) +
           rho / 2.0 * (F3 / F1) * ((F1 / rho) * (F1 / rho) + (F3 / F1) * (F3 / F1));
    return G;
}

void calculate_original(const vector<double>& F, double& rho, double& u, double& v, double& p, double& T) {
    double F1 = F[0], F2 = F[1], F3 = F[2], F4 = F[3];
    double A = F3 * F3 / (2.0 * F1) - F4;
    double B = GAMMA / (GAMMA - 1) * F1 * F2;
    double C = -(GAMMA + 1) / (2.0 * (GAMMA - 1)) * F1 * F1 * F1;
    rho = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    u = F1 / rho;
    v = F3 / F1;
    p = F2 - F1 * u;
    T = p / (rho * R);
}

void export_heatmap_data(const vector<double>& ksi_list, const vector<vector<double>>& y_field,
                         const vector<vector<double>>& Ma_field, const vector<vector<double>>& rho_field,
                         const vector<vector<double>>& p_field, const vector<vector<double>>& T_field, double rho1,
                         double p1, double T1);

int main() {
    vector<double> eta = linspace(0.0, max_eta, num_eta);
    vector<double> y_initial = linspace(0.0, H, num_eta);

    vector<double> ksi_list;
    vector<vector<double>> y_field;
    vector<vector<double>> rho_field, u_field, v_field, p_field, T_field, Ma_field;

    vector<double> rho(num_eta, rho1);
    vector<double> u(num_eta, Ma1 * sqrt(GAMMA * R * T1));
    vector<double> v(num_eta, 0.0);
    vector<double> p(num_eta, p1);
    vector<double> T(num_eta, T1);
    vector<double> Ma(num_eta, Ma1);

    vector<vector<double>> F(num_eta, vector<double>(4));
    for (int j = 0; j < num_eta; ++j) {
        F[j] = calculate_F(rho[j], u[j], v[j], p[j]);
    }

    ksi_list.push_back(0.0);
    y_field.push_back(y_initial);
    rho_field.push_back(rho);
    u_field.push_back(u);
    v_field.push_back(v);
    p_field.push_back(p);
    T_field.push_back(T);
    Ma_field.push_back(Ma);

    for (int iter = 0; iter < maxIteration; ++iter) {
        if (iter % 20 == 0) cout << "Iteration = " << iter << endl;

        double x_cur = ksi_list.back();
        if (x_cur >= L) break;

        double hx, P_eta_x0, angle0;
        height(x_cur, 0.0, hx, P_eta_x0, angle0);
        double dy = d_eta * hx;
        double max_tan = 0.0;
        for (int j = 0; j < num_eta; ++j) {
            double mu = asin(1.0 / Ma[j]);
            double tan1 = fabs(tan(angle0 + mu));
            double tan2 = fabs(tan(angle0 - mu));
            double m = max(tan1, tan2);
            if (m > max_tan) max_tan = m;
        }
        double dksi = C * dy / max_tan;
        vector<double> y_new(num_eta);
        for (int j = 0; j < num_eta; ++j) {
            y_new[j] = H - hx + j * dy;
        }
        y_field.push_back(y_new);

        // 预估步
        vector<vector<double>> G(num_eta, vector<double>(4));
        for (int j = 0; j < num_eta; ++j) {
            G[j] = calculate_G(F[j], rho[j]);
        }
        vector<vector<double>> P_F_ksi(num_eta, vector<double>(4, 0.0));
        for (int j = 0; j < num_eta - 1; ++j) {
            double P_eta_x;
            double h_tmp, angle_tmp;
            height(x_cur, eta[j], h_tmp, P_eta_x, angle_tmp);
            for (int k = 0; k < 4; ++k) {
                P_F_ksi[j][k] =
                    -(P_eta_x * (F[j + 1][k] - F[j][k]) / d_eta + 1.0 / hx * (G[j + 1][k] - G[j][k]) / d_eta);
            }
        }
        vector<vector<double>> SF(num_eta, vector<double>(4, 0.0));
        for (int j = 1; j < num_eta - 1; ++j) {
            double coeff = Cy * fabs(p[j + 1] - 2 * p[j] + p[j - 1]) / (p[j + 1] + 2 * p[j] + p[j - 1]);
            for (int k = 0; k < 4; ++k) {
                SF[j][k] = coeff * (F[j + 1][k] - 2 * F[j][k] + F[j - 1][k]);
            }
        }
        for (int k = 0; k < 4; ++k) {
            SF[0][k] = 2 * SF[1][k] - SF[2][k];
        }
        vector<vector<double>> F_pred(num_eta, vector<double>(4));
        for (int j = 0; j < num_eta; ++j) {
            for (int k = 0; k < 4; ++k) {
                F_pred[j][k] = F[j][k] + P_F_ksi[j][k] * dksi + SF[j][k];
            }
        }
        vector<double> rho_pred(num_eta), p_pred(num_eta);
        for (int j = 0; j < num_eta; ++j) {
            double rho_j, u_j, v_j, p_j, T_j;
            calculate_original(F_pred[j], rho_j, u_j, v_j, p_j, T_j);
            rho_pred[j] = rho_j;
            p_pred[j] = p_j;
        }

        // 校正步
        vector<vector<double>> G_pred(num_eta, vector<double>(4));
        for (int j = 0; j < num_eta; ++j) {
            G_pred[j] = calculate_G(F_pred[j], rho_pred[j]);
        }
        vector<vector<double>> P_F_ksi_pred(num_eta, vector<double>(4, 0.0));
        for (int j = 1; j < num_eta; ++j) {
            double P_eta_x;
            double h_tmp, angle_tmp;
            height(x_cur, eta[j], h_tmp, P_eta_x, angle_tmp);
            for (int k = 0; k < 4; ++k) {
                P_F_ksi_pred[j][k] = -(P_eta_x * (F_pred[j][k] - F[j - 1][k]) / d_eta +
                                       1.0 / hx * (G_pred[j][k] - G_pred[j - 1][k]) / d_eta);
            }
        }
        vector<vector<double>> P_av(num_eta, vector<double>(4));
        for (int j = 0; j < num_eta; ++j) {
            for (int k = 0; k < 4; ++k) {
                P_av[j][k] = 0.5 * (P_F_ksi[j][k] + P_F_ksi_pred[j][k]);
            }
        }
        vector<vector<double>> SF_pred(num_eta, vector<double>(4, 0.0));
        for (int j = 1; j < num_eta - 2; ++j) {
            double coeff = Cy * fabs(p_pred[j + 1] - 2 * p_pred[j] + p_pred[j - 1]) /
                           (p_pred[j + 1] + 2 * p_pred[j] + p_pred[j - 1]);
            for (int k = 0; k < 4; ++k) {
                SF_pred[j][k] = coeff * (F_pred[j + 1][k] - 2 * F_pred[j][k] + F_pred[j - 1][k]);
            }
        }
        for (int k = 0; k < 4; ++k) {
            SF_pred[num_eta - 2][k] = 2 * SF_pred[num_eta - 3][k] - SF_pred[num_eta - 4][k];
        }
        for (int j = 1; j < num_eta - 1; ++j) {
            for (int k = 0; k < 4; ++k) {
                F[j][k] = F[j][k] + P_av[j][k] * dksi + SF_pred[j][k];
            }
        }
        for (int k = 0; k < 4; ++k) {
            F[0][k] = 2 * F[1][k] - F[2][k];
        }

        for (int j = 0; j < num_eta; ++j) {
            double rho_j, u_j, v_j, p_j, T_j;
            calculate_original(F[j], rho_j, u_j, v_j, p_j, T_j);
            rho[j] = rho_j;
            u[j] = u_j;
            v[j] = v_j;
            p[j] = p_j;
            T[j] = T_j;
            Ma[j] = sqrt(u_j * u_j + v_j * v_j) / sqrt(GAMMA * R * T_j);
        }

        double angle_cur;
        double h_tmp, P_eta_x_tmp;
        height(x_cur, eta[0], h_tmp, P_eta_x_tmp, angle_cur);
        double phi = angle_cur + atan2(v[0], u[0]);
        AnalyticSolution sol = solution_analytic(Ma[0], p[0], T[0], rho[0], phi, angle_cur);
        Ma[0] = sol.Ma;
        p[0] = sol.p;
        T[0] = sol.T;
        rho[0] = sol.rho;
        u[0] = sol.u;
        v[0] = sol.v;
        F[0] = calculate_F(rho[0], u[0], v[0], p[0]);

        // 保存当前剖面
        ksi_list.push_back(x_cur + dksi);
        rho_field.push_back(rho);
        u_field.push_back(u);
        v_field.push_back(v);
        p_field.push_back(p);
        T_field.push_back(T);
        Ma_field.push_back(Ma);

        // 每隔10步输出一个单独的CSV文件
        int current_idx = ksi_list.size() - 1;  // 当前剖面的索引
        if (current_idx % 10 == 0 && current_idx != 0) {
            stringstream ss;
            ss << "../output/profile_" << current_idx << ".csv";
            ofstream fout(ss.str());
            fout << "x,y,Ma,rho,u,v,p,T\n";
            double x_out = ksi_list.back();
            const auto& y_out = y_field.back();
            const auto& Ma_out = Ma_field.back();
            const auto& rho_out = rho_field.back();
            const auto& u_out = u_field.back();
            const auto& v_out = v_field.back();
            const auto& p_out = p_field.back();
            const auto& T_out = T_field.back();
            for (int j = 0; j < num_eta; ++j) {
                fout << fixed << setprecision(3) << x_out << "," << y_out[j] << "," << Ma_out[j] << "," << rho_out[j]
                     << "," << u_out[j] << "," << v_out[j] << "," << p_out[j] << "," << T_out[j] << "\n";
            }
            fout.close();
            cout << "Saved profile " << current_idx << " to " << ss.str() << endl;
        }
    }

    cout << "Computation completed." << endl;
    // 导出热力图数据（所有剖面的所有网格点）
    // export_heatmap_data(ksi_list, y_field, Ma_field, rho_field, p_field, T_field, rho1, p1, T1);
    return 0;
}

void export_heatmap_data(const vector<double>& ksi_list, const vector<vector<double>>& y_field,
                         const vector<vector<double>>& Ma_field, const vector<vector<double>>& rho_field,
                         const vector<vector<double>>& p_field, const vector<vector<double>>& T_field, double rho1,
                         double p1, double T1) {
    // 打开输出文件（目录 ../output/ 需已存在）
    string filename = "../output/field_data.csv";
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    // 写入表头
    fout << "x,y,Ma,rho/rho1,p/p1,T/T1\n";

    // 设置输出格式
    fout << fixed << setprecision(6);

    int n_profiles = ksi_list.size();
    int n_eta = y_field[0].size();  // 网格点数（每个剖面的 y 点数目）

    // 遍历所有剖面和所有网格点
    for (int i = 0; i < n_profiles; ++i) {
        double x = ksi_list[i];
        const auto& y_vec = y_field[i];
        const auto& Ma_vec = Ma_field[i];
        const auto& rho_vec = rho_field[i];
        const auto& p_vec = p_field[i];
        const auto& T_vec = T_field[i];

        for (int j = 0; j < n_eta; ++j) {
            fout << x << "," << y_vec[j] << "," << Ma_vec[j] << "," << rho_vec[j] / rho1 << "," << p_vec[j] / p1 << ","
                 << T_vec[j] / T1 << "\n";
        }
    }

    fout.close();
    cout << "Heatmap data exported to " << filename << endl;
}