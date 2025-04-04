/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <eos/maths/power-of.hh>
#include <eos/nonleptonic-amplitudes/su3f-amplitudes.hh>
#include <eos/utils/options-impl.hh>

#include <map>

namespace eos
{
    using std::sqrt;

    NonleptonicAmplitudes<PToPP> *
    SU3FRepresentation<PToPP>::make(const Parameters & p, const Options & o)
    {
        return new SU3FRepresentation<PToPP>(p, o);
    }

    SU3FRepresentation<PToPP>::SU3FRepresentation(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        opt_q(o, options, "q"),
        opt_p1(o, options, "P1"),
        opt_p2(o, options, "P2"),
        opt_cp_conjugate(o, options, "cp-conjugate"),
        opt_B_bar(o, options, "B_bar"),
        theta_18(p["eta::theta_18"], *this),
        B(su3f::psd_b_triplet.find(opt_q.value())->second),
        P1{ {} },
        P2{ {} },
        Gfermi(p["WET::G_Fermi"], *this),
        re_AT3(p["nonleptonic::Re{AT3}@SU3F"], *this),
        im_AT3(p["nonleptonic::Im{AT3}@SU3F"], *this),
        re_CT3(p["nonleptonic::Re{CT3}@SU3F"], *this),
        im_CT3(p["nonleptonic::Im{CT3}@SU3F"], *this),
        re_AT6(p["nonleptonic::Re{AT6}@SU3F"], *this),
        im_AT6(p["nonleptonic::Im{AT6}@SU3F"], *this),
        re_CT6(p["nonleptonic::Re{CT6}@SU3F"], *this),
        im_CT6(p["nonleptonic::Im{CT6}@SU3F"], *this),
        re_AT15(p["nonleptonic::Re{AT15}@SU3F"], *this),
        im_AT15(p["nonleptonic::Im{AT15}@SU3F"], *this),
        re_CT15(p["nonleptonic::Re{CT15}@SU3F"], *this),
        im_CT15(p["nonleptonic::Im{CT15}@SU3F"], *this),
        re_BT3(p["nonleptonic::Re{BT3}@SU3F"], *this),
        im_BT3(p["nonleptonic::Im{BT3}@SU3F"], *this),
        re_BT6(p["nonleptonic::Re{BT6}@SU3F"], *this),
        im_BT6(p["nonleptonic::Im{BT6}@SU3F"], *this),
        re_BT15(p["nonleptonic::Re{BT15}@SU3F"], *this),
        im_BT15(p["nonleptonic::Im{BT15}@SU3F"], *this),
        re_DT3(p["nonleptonic::Re{DT3}@SU3F"], *this),
        im_DT3(p["nonleptonic::Im{DT3}@SU3F"], *this),
        re_AP3(p["nonleptonic::Re{AP3}@SU3F"], *this),
        im_AP3(p["nonleptonic::Im{AP3}@SU3F"], *this),
        re_CP3(p["nonleptonic::Re{CP3}@SU3F"], *this),
        im_CP3(p["nonleptonic::Im{CP3}@SU3F"], *this),
        re_AP6(p["nonleptonic::Re{AP6}@SU3F"], *this),
        im_AP6(p["nonleptonic::Im{AP6}@SU3F"], *this),
        re_CP6(p["nonleptonic::Re{CP6}@SU3F"], *this),
        im_CP6(p["nonleptonic::Im{CP6}@SU3F"], *this),
        re_AP15(p["nonleptonic::Re{AP15}@SU3F"], *this),
        im_AP15(p["nonleptonic::Im{AP15}@SU3F"], *this),
        re_CP15(p["nonleptonic::Re{CP15}@SU3F"], *this),
        im_CP15(p["nonleptonic::Im{CP15}@SU3F"], *this),
        re_BP3(p["nonleptonic::Re{BP3}@SU3F"], *this),
        im_BP3(p["nonleptonic::Im{BP3}@SU3F"], *this),
        re_BP6(p["nonleptonic::Re{BP6}@SU3F"], *this),
        im_BP6(p["nonleptonic::Im{BP6}@SU3F"], *this),
        re_BP15(p["nonleptonic::Re{BP15}@SU3F"], *this),
        im_BP15(p["nonleptonic::Im{BP15}@SU3F"], *this),
        re_DP3(p["nonleptonic::Re{DP3}@SU3F"], *this),
        im_DP3(p["nonleptonic::Im{DP3}@SU3F"], *this)
    {
        Context ctx("When constructing B->PP SU3 amplitudes");

        if (opt_cp_conjugate.value() != opt_B_bar.value())
        {
            lamdu = [this]() { return model->ckm_ub() * conj(model->ckm_ud()); };
            lamsu = [this]() { return model->ckm_ub() * conj(model->ckm_us()); };
            lamdt = [this]() { return model->ckm_tb() * conj(model->ckm_td()); };
            lamst = [this]() { return model->ckm_tb() * conj(model->ckm_ts()); };
        }
        else
        {
            lamdu = [this]() { return conj(model->ckm_ub()) * model->ckm_ud(); };
            lamsu = [this]() { return conj(model->ckm_ub()) * model->ckm_us(); };
            lamdt = [this]() { return conj(model->ckm_tb()) * model->ckm_td(); };
            lamst = [this]() { return conj(model->ckm_tb()) * model->ckm_ts(); };
        }

        H3bar = [this]()
        {
            return su3f::rank1{
                { 0.0, lamdu(), lamsu() }
            };
        };
        H3tilde = [this]()
        {
            return su3f::rank1{
                { 0.0, lamdt(), lamst() }
            };
        };

        H6bar = [this]()
        {
            su3f::rank3 result;
            result[0][1][0] = +lamdu();
            result[1][0][0] = -lamdu();
            result[1][2][2] = +lamdu();
            result[2][1][2] = -lamdu();
            result[0][2][0] = +lamsu();
            result[2][0][0] = -lamsu();
            result[2][1][1] = +lamsu(); // Corrected with respect to typo in [HTX:2021A]
            result[1][2][1] = -lamsu(); // Corrected with respect to typo in [HTX:2021A]
            return result;
        };

        H15bar = [this]()
        {
            su3f::rank3 result;
            result[0][1][0] = +3.0 * lamdu();
            result[1][0][0] = +3.0 * lamdu();
            result[1][1][1] = -2.0 * lamdu();
            result[1][2][2] = -lamdu();
            result[2][1][2] = -lamdu();
            result[0][2][0] = +3.0 * lamsu();
            result[2][0][0] = +3.0 * lamsu();
            result[2][2][2] = -2.0 * lamsu();
            result[2][1][1] = -lamsu(); // Corrected with respect to typo in [HTX:2021A]
            result[1][2][1] = -lamsu(); // Corrected with respect to typo in [HTX:2021A]
            return result;
        };

        H6tilde = [this]()
        {
            su3f::rank3 result;
            result[0][1][0] = +lamdt();
            result[1][0][0] = -lamdt();
            result[1][2][2] = +lamdt();
            result[2][1][2] = -lamdt();
            result[0][2][0] = +lamst();
            result[2][0][0] = -lamst();
            result[2][1][1] = +lamst(); // Corrected with respect to typo in [HTX:2021A]
            result[1][2][1] = -lamst(); // Corrected with respect to typo in [HTX:2021A]
            return result;
        };

        H15tilde = [this]()
        {
            su3f::rank3 result;
            result[0][1][0] = +3.0 * lamdt();
            result[1][0][0] = +3.0 * lamdt();
            result[1][1][1] = -2.0 * lamdt();
            result[1][2][2] = -lamdt();
            result[2][1][2] = -lamdt();
            result[0][2][0] = +3.0 * lamst();
            result[2][0][0] = +3.0 * lamst();
            result[2][2][2] = -2.0 * lamst();
            result[2][1][1] = -lamst(); // Corrected with respect to typo in [HTX:2021A]
            result[1][2][1] = -lamst(); // Corrected with respect to typo in [HTX:2021A]
            return result;
        };
    }

    const std::vector<OptionSpecification> SU3FRepresentation<PToPP>::options{
        Model::option_specification(),
        { "cp-conjugate",                                                                     { "true", "false" }, "false" },
        {        "B_bar",                                                                     { "true", "false" }, "false" },
        {            "q",                                                                       { "u", "d", "s" },      "" },
        {           "P1", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" },      "" },
        {           "P2", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" },      "" },
    };

    complex<double>
    SU3FRepresentation<PToPP>::tree_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const
    {
        complex<double>       T_ira = 0.0;
        const complex<double> AT3 = complex<double>(this->re_AT3(), this->im_AT3()), CT3 = complex<double>(this->re_CT3(), this->im_CT3()),
                              AT6 = complex<double>(this->re_AT6(), this->im_AT6()), CT6 = complex<double>(this->re_CT6(), this->im_CT6()),
                              AT15 = complex<double>(this->re_AT15(), this->im_AT15()), CT15 = complex<double>(this->re_CT15(), this->im_CT15()),
                              BT3 = complex<double>(this->re_BT3(), this->im_BT3()), BT6 = complex<double>(this->re_BT6(), this->im_BT6()),
                              BT15 = complex<double>(this->re_BT15(), this->im_BT15()), DT3 = complex<double>(this->re_DT3(), this->im_DT3());

        const auto & H3bar  = this->H3bar();
        const auto & H6bar  = this->H6bar();
        const auto & H15bar = this->H15bar();

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    T_ira += AT3 * B[i] * H3bar[i] * p1[j][k] * p2[k][j];
                    T_ira += CT3 * B[i] * p1[i][j] * p2[j][k] * H3bar[k];
                    T_ira += BT3 * B[i] * H3bar[i] * p1[k][k] * p2[j][j];
                    T_ira += DT3 * B[i] * p1[i][j] * H3bar[j] * p2[k][k];

                    for (unsigned l = 0; l < 3; l++)
                    {
                        T_ira += AT6 * B[i] * H6bar[i][j][k] * p1[l][j] * p2[k][l];
                        T_ira += CT6 * B[i] * p1[i][j] * H6bar[j][l][k] * p2[k][l];
                        T_ira += BT6 * B[i] * H6bar[i][j][k] * p1[k][j] * p2[l][l];

                        T_ira += AT15 * B[i] * H15bar[i][j][k] * p1[l][j] * p2[k][l];
                        T_ira += CT15 * B[i] * p1[i][j] * H15bar[j][k][l] * p2[l][k];
                        T_ira += BT15 * B[i] * H15bar[i][j][k] * p1[k][j] * p2[l][l];
                    }
                }
            }
        }

        return T_ira;
    }

    complex<double>
    SU3FRepresentation<PToPP>::penguin_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const
    {
        complex<double>       P_ira = 0.0;
        const complex<double> AP3 = complex<double>(this->re_AP3(), this->im_AP3()), CP3 = complex<double>(this->re_CP3(), this->im_CP3()),
                              AP6 = complex<double>(this->re_AP6(), this->im_AP6()), CP6 = complex<double>(this->re_CP6(), this->im_CP6()),
                              AP15 = complex<double>(this->re_AP15(), this->im_AP15()), CP15 = complex<double>(this->re_CP15(), this->im_CP15()),
                              BP3 = complex<double>(this->re_BP3(), this->im_BP3()), BP6 = complex<double>(this->re_BP6(), this->im_BP6()),
                              BP15 = complex<double>(this->re_BP15(), this->im_BP15()), DP3 = complex<double>(this->re_DP3(), this->im_DP3());

        const auto & H3tilde  = this->H3tilde();
        const auto & H6tilde  = this->H6tilde();
        const auto & H15tilde = this->H15tilde();

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    P_ira += AP3 * B[i] * H3tilde[i] * p1[j][k] * p2[k][j];
                    P_ira += CP3 * B[i] * p1[i][j] * p2[j][k] * H3tilde[k];
                    P_ira += BP3 * B[i] * H3tilde[i] * p1[k][k] * p2[j][j];
                    P_ira += DP3 * B[i] * p1[i][j] * H3tilde[j] * p2[k][k];

                    for (unsigned l = 0; l < 3; l++)
                    {
                        P_ira += AP6 * B[i] * H6tilde[i][j][k] * p1[l][j] * p2[k][l];
                        P_ira += CP6 * B[i] * p1[i][j] * H6tilde[j][l][k] * p2[k][l];
                        P_ira += BP6 * B[i] * H6tilde[i][j][k] * p1[k][j] * p2[l][l];

                        P_ira += AP15 * B[i] * H15tilde[i][j][k] * p1[l][j] * p2[k][l];
                        P_ira += CP15 * B[i] * p1[i][j] * H15tilde[j][k][l] * p2[l][k];
                        P_ira += BP15 * B[i] * H15tilde[i][j][k] * p1[k][j] * p2[l][l];
                    }
                }
            }
        }

        return P_ira;
    }

    complex<double>
    SU3FRepresentation<PToPP>::ordered_amplitude() const
    {
        this->update();

        if (opt_B_bar.value())
        {
            su3f::transpose(P1);
            su3f::transpose(P2);
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * (this->tree_amplitude(P1, P2) + this->penguin_amplitude(P1, P2));
    }

    complex<double>
    SU3FRepresentation<PToPP>::inverse_amplitude() const
    {
        this->update();

        if (opt_B_bar.value())
        {
            su3f::transpose(P1);
            su3f::transpose(P2);
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * (this->tree_amplitude(P2, P1) + this->penguin_amplitude(P2, P1));
    }

    complex<double>
    SU3FRepresentation<PToPP>::penguin_correction() const
    {
        this->update();

        auto penguin = (this->penguin_amplitude(P1, P2) + this->penguin_amplitude(P2, P1)) / lamdt();
        auto tree    = (this->tree_amplitude(P1, P2) + this->tree_amplitude(P2, P1)) / lamdu();

        return -penguin / (tree - penguin);
    }
} // namespace eos
