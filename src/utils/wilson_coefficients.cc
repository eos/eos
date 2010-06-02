/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/wilson_coefficients.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <vector>

namespace wf
{
    void calculate_wilson_coefficients(const double & mu, Parameters & parameters)
    {
        // cf. [CMM1997], Eq. (27), p. 10
        const static std::vector<double> a = { 6.0/23.0, -12.0/23.0, 0.4086, -0.4230, -0.8994, 0.1456 };
        // cf. [CMM1997], Eq. (28), p. 10
        const static std::vector<std::vector<double>> Ahat =
        {
            {  1.0,       -1.0,        0.0,     0.0,     0.0,     0.0    },
            {  2.0/3.0,    1.0/3.0,    0.0,     0.0,     0.0,     0.0    },
            {  2.0/63.0,  -1.0/27.0,  -0.0659,  0.0595, -0.0218,  0.0335 },
            {  1.0/21.0,   1.0/9.0,    0.0237, -0.0173, -0.1336, -0.0316 },
            { -1.0/126.0,  1.0/108.0,  0.0094, -0.0100,  0.0010, -0.0017 },
            { -1.0/84.0,  -1.0/36.0,   0.0108,  0.0163,  0.0103,  0.0023 },
        };
        // cf. [CMM1997], Eq. (29), p. 10
        const static std::vector<std::vector<double>> Bhat =
        {
            {  5.9606,  1.0951,  0.0,     0.0,     0.0,     0.0    },
            {  1.9737, -1.3650,  0.0,     0.0,     0.0,     0.0    },
            { -0.5409,  1.6332,  1.6406, -1.6702, -0.2576, -0.2250 },
            {  2.2203,  2.0265, -4.1830, -0.7135, -1.8215,  0.7996 },
            {  0.0400, -0.1860, -0.1669,  0.1887,  0.0201,  0.0304 },
            { -0.2614, -0.1918,  0.4197,  0.0295,  0.1474, -0.0640 },
        };
        // cf. [CMM1997], Eq. (30), p. 10
        const static std::vector<std::vector<double>> Bhatprime =
        {
            {  2.0394,  5.9049,  0.0,     0.0,     0.0,     0.0    },
            {  1.3596, -1.9683,  0.0,     0.0,     0.0,     0.0    },
            {  0.0647,  0.2187, -0.2979, -0.6218,  0.1880, -0.1318 },
            {  0.0971, -0.6561,  0.1071,  0.1806,  1.1520,  0.1242 },
            { -0.0162, -0.0547,  0.0423,  0.1041, -0.0085,  0.0067 },
            { -0.0243,  0.1640,  0.0489, -0.1700, -0.0889, -0.0091 },
        };
        // cf. [CMM1997], Eq. (31), p. 10
        const static std::vector<std::vector<double>> BhatE =
        {
            {  0.0,  0.0,  0.0,     0.0,     0.0,     0.0    },
            {  0.0,  0.0,  0.0,     0.0,     0.0,     0.0    },
            {  0.0,  0.0, -0.1933,  0.1579,  0.1428, -0.1074 },
            {  0.0,  0.0,  0.0695, -0.0459,  0.8752,  0.1012 },
            {  0.0,  0.0,  0.0274, -0.0264, -0.0064,  0.0055 },
            {  0.0,  0.0,  0.0317,  0.0432, -0.0675, -0.0074 },
        };

        double m_t = parameters["mass::t"];
        double m_W = parameters["mass::W"];

        double x = std::pow(m_t / m_W, 2.0);
        // cf. [CMM1997], Eq. 
        double E = x * (18.0 - 11.0 * x - x * x) / (12.0 * std::pow(1.0 - x, 3.0))
            + x * x * (15.0 - 16.0 * x + 4.0 * x * x) / (6.0 * std::pow(1.0 - x, 4.0)) * std::log(x)
            - 2.0 / 3.0 * std::log(x);

        double alpha_s = QCD::alpha_s(mu);
        double eta = QCD::alpha_s(m_W) / alpha_s;
        double nlo = alpha_s / (4.0 * M_PI);

        std::vector<double> result = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        for (unsigned i(0) ; i < result.size() ; ++i)
        {
            // cf. [CMM1997], Eq. (25), p. 10
            for (unsigned j(0) ; j < a.size() ; ++j)
            {
                result[i] += (Ahat[i][j] + nlo * (Bhat[i][j] + eta * (Bhatprime[i][j] + BhatE[i][j] * E))) * std::pow(eta, a[j]);
            }

            parameters["c" + std::string(1, char('1' + i))] = result[i];
        }
    }
}
