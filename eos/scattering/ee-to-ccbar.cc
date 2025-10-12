/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <vector>
#include <array>
#include <memory>

namespace eos
{
    using namespace std::literals::string_literals;

    template <>
    struct Implementation<EEToCCBar>
    {
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_eff;
        UsedParameter m_D0;
        UsedParameter m_Dp;

        bool assume_isospin;

        static const inline std::array<std::string, 3> resonance_names =
        {
            "J/psi", "psi(2S)", "psi(3770)"
        };

        enum Resonances
        {
            Jpsi = 0, psi2S, psi3770
        };

        static const inline std::array<std::string, 6> channel_names =
        {
            "e^+e^-", "eff(Jpsi)", "eff(2S)", "D^0Dbar^0", "D^+D^-", "eff(3770)"
        };

        enum Channels
        {
            ee = 0, effJpsi, eff2S, D0Dbar0, DpDm, eff3770
        };

        // Resonance masses
        std::array<UsedParameter, EEToCCBar::nresonances> m;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nresonances> g0;

        // Channels barrier factors scales
        std::array<UsedParameter, EEToCCBar::nchannels> q;

        // Non resonant contributions to the K matrix
        std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nchannels> bkgcst;

        // R_uds
        UsedParameter Rconstant;

        // Normalization of the exclusive channels
        UsedParameter exclusive_norm;

        std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>> K;

        using IntermediateResult = EEToCCBar::IntermediateResult;
        IntermediateResult _intermediate_result;

        template <typename T, T... indices>
        auto _resonance_masses(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            static_assert(sizeof...(indices) <= resonance_names.size(), "The number of requested resonances is larger than the number of resonance names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                UsedParameter(p["mass::" + resonance_names[indices]], u)...
            }};
        }

        template <typename T, T... indices>
        auto _channel_effective_momentum(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            static_assert(sizeof...(indices) <= channel_names.size(), "The number of requested channels is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                // TODO: The q0 widths can be made channel dependent.
                // But in order to fit them all with one single parameter, we use a common parameter "ee->ccbar::q_0".
                // Once EOS provides a way of aliasing parameters, we can use the following line
                // UsedParameter(p["ee->ccbar::q_0(" + channel_names[indices] + ")"], u)...
                // Here we use a trick to allow C++ to expand the indices: the void followed by a coma is ignored,
                // and this line only instentiate an array of `UsedParameter(p["ee->ccbar::q_0"], u)`
                (static_cast<void>(indices), UsedParameter(p["ee->ccbar::q_0"], u))...
            }};
        }

        long unsigned _filter_channel_index(Channels channel)
        {
            if (assume_isospin)
            {
                switch (channel)
                {
                    case DpDm:
                        return D0Dbar0;

                    default:
                        return channel;
                }
            }
            else
            {
                switch (channel)
                {
                    default:
                        return channel;
                }
            }
        }

        template <typename T, T... column_indices>
        auto _g0_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            static_assert(sizeof...(column_indices) <= channel_names.size(), "The number of requested channels is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::g0(" + resonance_names[row_index] + "," + channel_names[_filter_channel_index(Channels(column_indices))] + ")"], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _g0_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            static_assert(sizeof...(row_indices) <= resonance_names.size(), "The number of requested resonances is larger than the number of resonance names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _g0_row(p, u, row_indices, column_seq)...
            }};
        }

        template <typename T>
        std::string _channel_name_tuple(const T & a, const T & b)
        {
            if (b > a)
                return "(" + channel_names[_filter_channel_index(Channels(a))] + "," + channel_names[_filter_channel_index(Channels(b))] + ")";

            return "(" + channel_names[_filter_channel_index(Channels(b))] + "," + channel_names[_filter_channel_index(Channels(a))] + ")";
        }

        template <typename T, T... column_indices>
        auto _c_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            static_assert(sizeof...(column_indices) <= channel_names.size(), "The number of requested channels (column) is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::c" + _channel_name_tuple(row_index, column_indices)], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _c_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            static_assert(sizeof...(row_indices) <= channel_names.size(), "The number of requested channels (row) is larger than the number of channel names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _c_row(p, u, row_indices, column_seq)...
            }};
        }

        template <typename T, T... row_indices>
        auto _get_g0_column(T column_index, std::integer_sequence<T, row_indices...>)
            -> std::array<Parameter, sizeof...(row_indices)>
        {
            return std::array<Parameter, sizeof...(row_indices)>
            {{
                g0[row_indices][column_index]...
            }};
        }

        template <typename T, T... column_indices>
        auto _get_bkgcst_row(T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<Parameter, sizeof...(column_indices)>
        {
            return std::array<Parameter, sizeof...(column_indices)>
            {{
                bkgcst[row_index][column_indices]...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _get_bkgcst_matrix(std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            return std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _get_bkgcst_row(row_indices, column_seq)...
            }};
        }


        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            hbar(p["QM::hbar"], u),
            alpha_em(p["QED::alpha_e(0)"], u),
            m_e(p["mass::e"], u),
            m_eff(p["ee->ccbar::effective_mass"], u),
            m_D0(p["mass::D^0"], u),
            m_Dp(p["mass::D^+"], u),
            assume_isospin(destringify<bool>(o.get("assume-isospin"_ok, "false"))),
            m(_resonance_masses(p, u, std::make_index_sequence<EEToCCBar::nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<EEToCCBar::nresonances>(), std::make_index_sequence<EEToCCBar::nchannels>())),
            q(_channel_effective_momentum(p, u, std::make_index_sequence<EEToCCBar::nchannels>())),
            bkgcst(_c_matrix(p, u, std::make_index_sequence<EEToCCBar::nchannels>(), std::make_index_sequence<EEToCCBar::nchannels>())),
            Rconstant(p["ee->ccbar::Rconstant"], u),
            exclusive_norm(p["ee->ccbar::exclusive_norm"], u)
        {
            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>::Resonance>, EEToCCBar::nresonances> resonance_array;
            for (unsigned i = 0 ; i < EEToCCBar::nresonances ; i++)
            {
                resonance_array[i] = std::make_shared<CharmoniumResonance<EEToCCBar::nchannels, EEToCCBar::nresonances>>(resonance_names[i], m[i]);
            }

            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>::Channel>, EEToCCBar::nchannels> channel_array;
            for (unsigned i = 0 ; i < EEToCCBar::nchannels ; i++)
            {
                switch (Channels(i))
                {
                    case ee:
                        channel_array[i] = std::make_shared<EEChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_e, m_e, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case effJpsi:
                    case eff2S:
                    case eff3770:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_eff, m_eff, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case D0Dbar0:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_D0, m_D0, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DpDm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dp, m_Dp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of known channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>> (
                new KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>(
                    channel_array,
                    resonance_array,
                    _get_bkgcst_matrix(std::make_index_sequence<EEToCCBar::nchannels>(), std::make_index_sequence<EEToCCBar::nchannels>()),
                    "e^+e^-->ccbar"
                    )
                );
        }

        const IntermediateResult * prepare(const complex<double> & E)
        {
            // Amplitude on the first RS
            _intermediate_result.tmatrix_row_0 = K->tmatrix_row(0, E * E);
            // Amplitude on the second RS
            _intermediate_result.tmatrix2_row_0 = K->tmatrix_row(0, E * E, true);

            _intermediate_result.E = E;
            _intermediate_result.s = E * E;

            return &_intermediate_result;
        }


        double rho(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            return std::real(K->_channels[channel]->rho(intermediate_result->s));
        }

        complex<double> chew_mandelstam(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            return K->_channels[channel]->chew_mandelstam(intermediate_result->s);
        }

        complex<double> chew_mandelstam_II(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            complex<double> s = intermediate_result->s;
            auto get_channel = K->_channels[channel];

            const unsigned li = get_channel->_l_orbital;
            const double q0 = get_channel->_q0.evaluate();
            const complex<double> mi1_2 = power_of<2>(get_channel->_m1());
            const complex<double> mi2_2 = power_of<2>(get_channel->_m2());

            // Momentum of particles in their center-of-momentum frame
            const complex<double> q = 0.5 * sqrt(eos::lambda(s, mi1_2, mi2_2)) / sqrt(s);

            // Blatt-Weisskopf factors, cf eq. (50.26)
            const complex<double> Fi = kmatrix_utils::blatt_weisskopf_factor(li, q / q0);

            return get_channel->chew_mandelstam(s) +
                complex<double>(0.0, 2.0) * get_channel->rho(s) * power_of<2>(pow(q / q0, li) * Fi);
        }


        inline double sigma_eetomumu_leading_order(const double & E)
        {
            // Conversion factor between GeV^2 and nb
            static const double speedoflight = 299792458.; // Exact value
            const double GeVtonb = 10.0 * power_of<2>(1.0e18 * hbar * speedoflight);

            return GeVtonb * 4.0 * M_PI * alpha_em * alpha_em / (3.0 * E * E);
        }

        // amplitude of ee -> channel
        complex<double> T_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[channel];

            return T1f;
        }

        complex<double> T_II_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Get T-matrix[ee, channel] on the second RS
            const complex<double> T1f = intermediate_result->tmatrix2_row_0[channel];

            return T1f;
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Conversion factor between GeV^2 and nb
            static const double speedoflight = 299792458.; // Exact value
            const double GeVtonb = 10.0 * power_of<2>(1.0e18 * hbar * speedoflight);

            // Channel properties
            const double Nf = 2.0 * K->_channels[channel]->_l_orbital + 1.0;
            const double rhof = std::abs(K->_channels[channel]->rho(intermediate_result->s));

            return GeVtonb / std::abs(intermediate_result->s) * Nf * rhof * norm(T_eetochannel(intermediate_result, channel));
        }

        // K matrix widths, they are not expected to match the experimental ones
        double res_partial_width(const Resonances & resonance, const Channels & channel)
        {
            return K->partial_width(resonance, channel);
        }

        double res_total_width(const Resonances & resonance)
        {
            return K->width(resonance);
        }

        // Ratios
        double R(const IntermediateResult * intermediate_result)
        {
            double total_sigma = 0.0;

            for (unsigned i = 1 ; i < EEToCCBar::nchannels ; i++)
            {
                total_sigma += sigma_eetochannel(intermediate_result, Channels(i));
            }

            return total_sigma / sigma_eetomumu_leading_order(abs(intermediate_result->E)) + Rconstant; // Add constant term
        }

        // Spectral function
        double spectral_function(const double E, const Resonances & resonance)
        {
            return K->spectral_function(resonance, E * E);
        }
    };

    EEToCCBar::EEToCCBar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBar>(new Implementation<EEToCCBar>(parameters, options, *this))
    {
    }

    EEToCCBar::~EEToCCBar()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<EEToCCBar>::options
    {
        {"assume-isospin"_ok, { "true"s, "false"s }, "false"s},
    };

    const EEToCCBar::IntermediateResult *
    EEToCCBar::prepare(const double & E) const
    {
        return _imp->prepare(E);
    }

    const EEToCCBar::IntermediateResult *
    EEToCCBar::prepare_complex(const double & re_E, const double & im_E) const
    {
        return _imp->prepare(complex<double>(re_E, im_E));
    }

    using Resonances = Implementation<eos::EEToCCBar>::Resonances;
    using Channels = Implementation<eos::EEToCCBar>::Channels;

    double
    EEToCCBar::Jpsi_ee_width() const
    {
        return _imp->res_partial_width(Resonances::Jpsi, Channels::ee);
    }

    double
    EEToCCBar::Jpsi_eff_width() const
    {
        return _imp->res_partial_width(Resonances::Jpsi, Channels::eff2S);
    }

    double
    EEToCCBar::Jpsi_total_width() const
    {
        return _imp->res_total_width(Resonances::Jpsi);
    }

    double
    EEToCCBar::psi2S_ee_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::ee);
    }

    double
    EEToCCBar::psi2S_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::eff2S);
    }

    double
    EEToCCBar::psi2S_total_width() const
    {
        return _imp->res_total_width(Resonances::psi2S);
    }

    double
    EEToCCBar::psi3770_total_width() const
    {
        return _imp->res_total_width(Resonances::psi3770);
    }

    double
    EEToCCBar::psi3770_D0Dbar0_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::D0Dbar0);
    }

    double
    EEToCCBar::psi3770_DpDm_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::DpDm);
    }

    double
    EEToCCBar::psi3770_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::eff3770);
    }

    double
    EEToCCBar::sigma_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::ee);
    }

    double
    EEToCCBar::sigma_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::effJpsi)
            + _imp->sigma_eetochannel(ir, Channels::eff2S)
            + _imp->sigma_eetochannel(ir, Channels::eff3770)
            );
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::D0Dbar0);
    }

    double
    EEToCCBar::sigma_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::DpDm);
    }

    double
    EEToCCBar::rho_ee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::ee);
    }

    double
    EEToCCBar::rho_eff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::eff3770);
    }

    double
    EEToCCBar::rho_D0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::D0Dbar0);
    }

    double
    EEToCCBar::rho_DpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::DpDm);
    }

    double
    EEToCCBar::re_chew_mandelstam_ee(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::ee));
    }

    double
    EEToCCBar::im_chew_mandelstam_ee(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::ee));
    }

    double
    EEToCCBar::re_chew_mandelstam_eff(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::eff3770));
    }

    double
    EEToCCBar::im_chew_mandelstam_eff(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::eff3770));
    }

    double
    EEToCCBar::re_chew_mandelstam_D0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::im_chew_mandelstam_D0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::re_chew_mandelstam_DpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::DpDm));
    }

    double
    EEToCCBar::im_chew_mandelstam_DpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::DpDm));
    }
    double
    EEToCCBar::re_chew_mandelstam_II_ee(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::ee));
    }

    double
    EEToCCBar::im_chew_mandelstam_II_ee(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::ee));
    }

    double
    EEToCCBar::re_chew_mandelstam_II_eff(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::eff3770));
    }

    double
    EEToCCBar::im_chew_mandelstam_II_eff(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::eff3770));
    }

    double
    EEToCCBar::re_chew_mandelstam_II_D0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::im_chew_mandelstam_II_D0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::re_chew_mandelstam_II_DpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::DpDm));
    }

    double
    EEToCCBar::im_chew_mandelstam_II_DpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::DpDm));
    }


    double
    EEToCCBar::re_T_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBar::im_T_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBar::re_T_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBar::im_T_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBar::re_T_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBar::im_T_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBar::re_T_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::im_T_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::re_T_II_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBar::im_T_II_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBar::re_T_II_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBar::im_T_II_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBar::re_T_II_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBar::im_T_II_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBar::re_T_II_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::im_T_II_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBar::psi3770_spectral_function(const double & E) const
    {
        return _imp->spectral_function(E, Resonances::psi3770);
    }


    double
    EEToCCBar::R(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->R(ir);
    }

    const std::set<ReferenceName>
    EEToCCBar::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    EEToCCBar::begin_options()
    {
        return Implementation<EEToCCBar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    EEToCCBar::end_options()
    {
        return Implementation<EEToCCBar>::options.cend();
    }
}
