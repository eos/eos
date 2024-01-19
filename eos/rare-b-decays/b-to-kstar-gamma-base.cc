/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/rare-b-decays/b-to-kstar-gamma-base.hh>
#include <eos/utils/destringify.hh>

namespace eos
{
    BToKstarGamma::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors", "BSZ2015"), p)),
        hbar(p["QM::hbar"], *this),
        mu(p["sb::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        m_B(p["mass::B_" + o.get("q", "d")], *this),
        m_Kstar(p["mass::K_d^*"], *this),
        l(o, options, "l"),
        m_l(p["mass::" + l.str()], *this),
        cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
    {
        std::string spectator_quark = o.get("q", "d");
        if (spectator_quark.size() != 1)
            throw InternalError("Option q should only be one character!");

        q = spectator_quark[0];
        if (q == 'd')
        {
            e_q = -1.0 / 3.0;
        }
        else if (q == 'u')
        {
            e_q = 2.0 / 3.0;
        }
        else
        {
            throw InternalError("Unsupported spectator quark");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKstarGamma::AmplitudeGenerator::~AmplitudeGenerator() = default;

    const std::vector<OptionSpecification>
    BToKstarGamma::AmplitudeGenerator::options
    {
        Model::option_specification(),
        { "l", { "e", "mu" }, "mu" },
    };
}
