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
        q(o, options, "q"),
        tau(p["life_time::B_" + q.str()], *this),
        m_B(p["mass::B_" + q.str()], *this),
        m_Kstar(p["mass::K_d^*"], *this),
        l(o, options, "l"),
        m_l(p["mass::" + l.str()], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"),
        cp_conjugate(opt_cp_conjugate.value())
    {
        Context ctx("When constructing B->K^*gamma amplitudes");

        switch (q.value())
        {
            case QuarkFlavor::down:
                e_q = -1.0 / 3.0;
                break;

            case QuarkFlavor::up:
                e_q = 2.0 / 3.0;
                break;

            default:
                throw InternalError("Unexpected quark flavor: '" + q.str() + "'");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKstarGamma::AmplitudeGenerator::~AmplitudeGenerator() = default;

    const std::vector<OptionSpecification>
    BToKstarGamma::AmplitudeGenerator::options
    {
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        { "l", { "e", "mu" }, "mu" },
        { "q", { "d", "u" }, "d" },
        { "cp-conjugate", { "true", "false" }, "false" }
    };
}
