/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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
#include <eos/observable-impl.hh>
#include <eos/observable.hh>
#include <eos/signal-pdf-impl.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/instantiation_policy-impl.hh>

#include <test/test-pdfs.hh>

#include <algorithm>
#include <cmath>

namespace eos::test
{
    namespace
    {
        // integral of exp(-(v - mu)^2 / (2 sigma^2)) over [lo, hi]
        double
        gaussian_norm(const double & mu, const double & sigma, const double & lo, const double & hi)
        {
            const double s = sigma * M_SQRT2;

            return sigma * std::sqrt(M_PI / 2.0) * (std::erf((hi - mu) / s) - std::erf((lo - mu) / s));
        }

        // length of the overlap of [centre - width / 2, centre + width / 2] with [lo, hi]
        double
        box_norm(const double & centre, const double & width, const double & lo, const double & hi)
        {
            const double a = std::max(lo, centre - 0.5 * width);
            const double b = std::min(hi, centre + 0.5 * width);

            return std::max(0.0, b - a);
        }
    } // namespace

    // Unnormalized PDF = exp(-(x - mu)^2 / (2 sigma^2))
    class TestGaussian1DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter mu;

            UsedParameter sigma;

            TestGaussian1DPDF(const Parameters & p, const Options &) :
                mu(p["TestGaussian1D::mu"], *this),
                sigma(p["TestGaussian1D::sigma"], *this)
            {
            }

            double
            pdf(const double & x) const
            {
                return std::exp(-power_of<2>(x - mu()) / (2.0 * power_of<2>(sigma())));
            }

            double
            norm(const double & x_min, const double & x_max) const
            {
                return gaussian_norm(mu(), sigma(), x_min, x_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestGaussian1DPDF::options{};

    const std::set<ReferenceName> TestGaussian1DPDF::references{};

    const std::string TestGaussian1DPDF::description = "1D Gaussian PDF as a function of x; used for unit tests only.";

    // Unnormalized PDF = 1 inside [centre - width / 2, centre + width / 2], else 0
    class TestBox1DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter centre;

            UsedParameter width;

            TestBox1DPDF(const Parameters & p, const Options &) :
                centre(p["TestBox1D::centre"], *this),
                width(p["TestBox1D::width"], *this)
            {
            }

            double
            pdf(const double & x) const
            {
                return (std::abs(x - centre()) <= 0.5 * width()) ? 1.0 : 0.0;
            }

            double
            norm(const double & x_min, const double & x_max) const
            {
                return box_norm(centre(), width(), x_min, x_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestBox1DPDF::options{};

    const std::set<ReferenceName> TestBox1DPDF::references{};

    const std::string TestBox1DPDF::description = "1D top-hat PDF as a function of x; used for unit tests only.";

    // Unnormalized PDF = product of two uncorrelated 1D Gaussians
    class TestGaussian2DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter mu_x;

            UsedParameter sigma_x;

            UsedParameter mu_y;

            UsedParameter sigma_y;

            TestGaussian2DPDF(const Parameters & p, const Options &) :
                mu_x(p["TestGaussian2D::mu_x"], *this),
                sigma_x(p["TestGaussian2D::sigma_x"], *this),
                mu_y(p["TestGaussian2D::mu_y"], *this),
                sigma_y(p["TestGaussian2D::sigma_y"], *this)
            {
            }

            double
            pdf(const double & x, const double & y) const
            {
                return std::exp(-power_of<2>(x - mu_x()) / (2.0 * power_of<2>(sigma_x())) - power_of<2>(y - mu_y()) / (2.0 * power_of<2>(sigma_y())));
            }

            double
            norm(const double & x_min, const double & x_max, const double & y_min, const double & y_max) const
            {
                return gaussian_norm(mu_x(), sigma_x(), x_min, x_max) * gaussian_norm(mu_y(), sigma_y(), y_min, y_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestGaussian2DPDF::options{};

    const std::set<ReferenceName> TestGaussian2DPDF::references{};

    const std::string TestGaussian2DPDF::description = "2D PDF as a function of x, y; a separable product of two 1D Gaussians; used for unit tests only.";

    // Unnormalized PDF = product of two 1D top-hats
    class TestBox2DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter centre_x;

            UsedParameter width_x;

            UsedParameter centre_y;

            UsedParameter width_y;

            TestBox2DPDF(const Parameters & p, const Options &) :
                centre_x(p["TestBox2D::centre_x"], *this),
                width_x(p["TestBox2D::width_x"], *this),
                centre_y(p["TestBox2D::centre_y"], *this),
                width_y(p["TestBox2D::width_y"], *this)
            {
            }

            double
            pdf(const double & x, const double & y) const
            {
                return ((std::abs(x - centre_x()) <= 0.5 * width_x()) && (std::abs(y - centre_y()) <= 0.5 * width_y())) ? 1.0 : 0.0;
            }

            double
            norm(const double & x_min, const double & x_max, const double & y_min, const double & y_max) const
            {
                return box_norm(centre_x(), width_x(), x_min, x_max) * box_norm(centre_y(), width_y(), y_min, y_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestBox2DPDF::options{};

    const std::set<ReferenceName> TestBox2DPDF::references{};

    const std::string TestBox2DPDF::description = "2D PDF as a function of x, y; a separable product of two 1D top-hats; used for unit tests only.";

    // Unnormalized PDF = exp(-x^2 / (2 sigma^2)); a Gaussian resolution over an offset variable
    class TestGaussianResolution1DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter sigma;

            TestGaussianResolution1DPDF(const Parameters & p, const Options &) :
                sigma(p["TestGaussianResolution1D::sigma"], *this)
            {
            }

            double
            pdf(const double & x) const
            {
                return std::exp(-power_of<2>(x) / (2.0 * power_of<2>(sigma())));
            }

            double
            norm(const double & x_min, const double & x_max) const
            {
                return gaussian_norm(0.0, sigma(), x_min, x_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestGaussianResolution1DPDF::options{};

    const std::set<ReferenceName> TestGaussianResolution1DPDF::references{};

    const std::string TestGaussianResolution1DPDF::description = "1D Gaussian resolution PDF as a function of the offset x; used for unit tests only.";

    // Unnormalized PDF = 1 inside [-width / 2, width / 2], else 0; a top-hat resolution over an offset variable
    class TestBoxResolution1DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter width;

            TestBoxResolution1DPDF(const Parameters & p, const Options &) :
                width(p["TestBoxResolution1D::width"], *this)
            {
            }

            double
            pdf(const double & x) const
            {
                return (std::abs(x) <= 0.5 * width()) ? 1.0 : 0.0;
            }

            double
            norm(const double & x_min, const double & x_max) const
            {
                return box_norm(0.0, width(), x_min, x_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestBoxResolution1DPDF::options{};

    const std::set<ReferenceName> TestBoxResolution1DPDF::references{};

    const std::string TestBoxResolution1DPDF::description = "1D top-hat resolution PDF as a function of the offset x; used for unit tests only.";

    // Unnormalized PDF = exp(-x^2 / (2 sigma_x^2) - y^2 / (2 sigma_y^2)); a separable 2D Gaussian resolution over the offsets
    class TestGaussianResolution2DPDF : public ParameterUser
    {
        public:
            static const std::vector<OptionSpecification> options;

            static const std::set<ReferenceName> references;

            static const std::string description;

            UsedParameter sigma_x;

            UsedParameter sigma_y;

            TestGaussianResolution2DPDF(const Parameters & p, const Options &) :
                sigma_x(p["TestGaussianResolution2D::sigma_x"], *this),
                sigma_y(p["TestGaussianResolution2D::sigma_y"], *this)
            {
            }

            double
            pdf(const double & x, const double & y) const
            {
                return std::exp(-power_of<2>(x) / (2.0 * power_of<2>(sigma_x())) - power_of<2>(y) / (2.0 * power_of<2>(sigma_y())));
            }

            double
            norm(const double & x_min, const double & x_max, const double & y_min, const double & y_max) const
            {
                return gaussian_norm(0.0, sigma_x(), x_min, x_max) * gaussian_norm(0.0, sigma_y(), y_min, y_max);
            }

            static ObservableEntry::OptionIterator
            begin_options()
            {
                return options.begin();
            }

            static ObservableEntry::OptionIterator
            end_options()
            {
                return options.end();
            }
    };

    const std::vector<OptionSpecification> TestGaussianResolution2DPDF::options{};

    const std::set<ReferenceName> TestGaussianResolution2DPDF::references{};

    const std::string TestGaussianResolution2DPDF::description = "2D Gaussian resolution PDF as a function of the offsets x, y; used for unit tests only.";

    void
    register_test_pdfs()
    {
        static bool registered = false;
        if (registered)
        {
            return;
        }
        registered = true;

        using std::literals::string_literals::operator""s;

        Parameters::declare("TestGaussian1D::mu", R"(\mu)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestGaussian1D::sigma", R"(\sigma)", Unit::None(), 1.0, 1.0e-6, 100.0);
        Parameters::declare("TestBox1D::centre", R"(c)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestBox1D::width", R"(w)", Unit::None(), 2.0, 1.0e-6, 100.0);
        Parameters::declare("TestGaussian2D::mu_x", R"(\mu_x)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestGaussian2D::sigma_x", R"(\sigma_x)", Unit::None(), 1.0, 1.0e-6, 100.0);
        Parameters::declare("TestGaussian2D::mu_y", R"(\mu_y)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestGaussian2D::sigma_y", R"(\sigma_y)", Unit::None(), 1.0, 1.0e-6, 100.0);
        Parameters::declare("TestBox2D::centre_x", R"(c_x)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestBox2D::width_x", R"(w_x)", Unit::None(), 2.0, 1.0e-6, 100.0);
        Parameters::declare("TestBox2D::centre_y", R"(c_y)", Unit::None(), 0.0, -100.0, 100.0);
        Parameters::declare("TestBox2D::width_y", R"(w_y)", Unit::None(), 2.0, 1.0e-6, 100.0);
        Parameters::declare("TestGaussianResolution1D::sigma", R"(\sigma_{\rm res})", Unit::None(), 1.0, 1.0e-6, 100.0);
        Parameters::declare("TestBoxResolution1D::width", R"(w_{\rm res})", Unit::None(), 2.0, 1.0e-6, 100.0);
        Parameters::declare("TestGaussianResolution2D::sigma_x", R"(\sigma_{{\rm res},x})", Unit::None(), 1.0, 1.0e-6, 100.0);
        Parameters::declare("TestGaussianResolution2D::sigma_y", R"(\sigma_{{\rm res},y})", Unit::None(), 1.0, 1.0e-6, 100.0);

        // TestGaussian1D::P(x)
        {
            make_observable("TestGaussian1D::UnnormalizedPDF(x)", Unit::None(), &TestGaussian1DPDF::pdf, std::make_tuple("x"));
            make_observable("TestGaussian1D::NormalizationPDF(x)", Unit::None(), &TestGaussian1DPDF::norm, std::make_tuple("x_min", "x_max"));

            auto signal_pdf = make_signal_pdf("TestGaussian1D::P(x)",
                                              "PDF for testing purpose only: 1D Gaussian PDF as a function of x.",
                                              Options{},
                                              "TestGaussian1D::UnnormalizedPDF(x)",
                                              std::make_tuple("x"s),
                                              "TestGaussian1D::NormalizationPDF(x)",
                                              std::make_tuple("x_min"s, "x_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestBox1D::P(x)
        {
            make_observable("TestBox1D::UnnormalizedPDF(x)", Unit::None(), &TestBox1DPDF::pdf, std::make_tuple("x"));
            make_observable("TestBox1D::NormalizationPDF(x)", Unit::None(), &TestBox1DPDF::norm, std::make_tuple("x_min", "x_max"));

            auto signal_pdf = make_signal_pdf("TestBox1D::P(x)",
                                              "PDF for testing purpose only: 1D top-hat PDF as a function of x.",
                                              Options{},
                                              "TestBox1D::UnnormalizedPDF(x)",
                                              std::make_tuple("x"s),
                                              "TestBox1D::NormalizationPDF(x)",
                                              std::make_tuple("x_min"s, "x_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestGaussian2D::P(x,y)
        {
            make_observable("TestGaussian2D::UnnormalizedPDF(x,y)", Unit::None(), &TestGaussian2DPDF::pdf, std::make_tuple("x", "y"));
            make_observable("TestGaussian2D::NormalizationPDF(x,y)", Unit::None(), &TestGaussian2DPDF::norm, std::make_tuple("x_min", "x_max", "y_min", "y_max"));

            auto signal_pdf = make_signal_pdf("TestGaussian2D::P(x,y)",
                                              "PDF for testing purpose only: 2D PDF as a function of x, y, separable into two 1D Gaussians.",
                                              Options{},
                                              "TestGaussian2D::UnnormalizedPDF(x,y)",
                                              std::make_tuple("x"s, "y"s),
                                              "TestGaussian2D::NormalizationPDF(x,y)",
                                              std::make_tuple("x_min"s, "x_max"s, "y_min"s, "y_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestBox2D::P(x,y)
        {
            make_observable("TestBox2D::UnnormalizedPDF(x,y)", Unit::None(), &TestBox2DPDF::pdf, std::make_tuple("x", "y"));
            make_observable("TestBox2D::NormalizationPDF(x,y)", Unit::None(), &TestBox2DPDF::norm, std::make_tuple("x_min", "x_max", "y_min", "y_max"));

            auto signal_pdf = make_signal_pdf("TestBox2D::P(x,y)",
                                              "PDF for testing purpose only: 2D PDF as a function of x, y, separable into two 1D top-hats.",
                                              Options{},
                                              "TestBox2D::UnnormalizedPDF(x,y)",
                                              std::make_tuple("x"s, "y"s),
                                              "TestBox2D::NormalizationPDF(x,y)",
                                              std::make_tuple("x_min"s, "x_max"s, "y_min"s, "y_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestGaussianResolution1D::P(x)
        {
            make_observable("TestGaussianResolution1D::UnnormalizedPDF(x)", Unit::None(), &TestGaussianResolution1DPDF::pdf, std::make_tuple("x"));
            make_observable("TestGaussianResolution1D::NormalizationPDF(x)", Unit::None(), &TestGaussianResolution1DPDF::norm, std::make_tuple("x_min", "x_max"));

            auto signal_pdf = make_signal_pdf("TestGaussianResolution1D::P(x)",
                                              "PDF for testing purpose only: 1D Gaussian resolution as a function of the offset x.",
                                              Options{},
                                              "TestGaussianResolution1D::UnnormalizedPDF(x)",
                                              std::make_tuple("x"s),
                                              "TestGaussianResolution1D::NormalizationPDF(x)",
                                              std::make_tuple("x_min"s, "x_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestBoxResolution1D::P(x)
        {
            make_observable("TestBoxResolution1D::UnnormalizedPDF(x)", Unit::None(), &TestBoxResolution1DPDF::pdf, std::make_tuple("x"));
            make_observable("TestBoxResolution1D::NormalizationPDF(x)", Unit::None(), &TestBoxResolution1DPDF::norm, std::make_tuple("x_min", "x_max"));

            auto signal_pdf = make_signal_pdf("TestBoxResolution1D::P(x)",
                                              "PDF for testing purpose only: 1D top-hat resolution as a function of the offset x.",
                                              Options{},
                                              "TestBoxResolution1D::UnnormalizedPDF(x)",
                                              std::make_tuple("x"s),
                                              "TestBoxResolution1D::NormalizationPDF(x)",
                                              std::make_tuple("x_min"s, "x_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }

        // TestGaussianResolution2D::P(x,y)
        {
            make_observable("TestGaussianResolution2D::UnnormalizedPDF(x,y)", Unit::None(), &TestGaussianResolution2DPDF::pdf, std::make_tuple("x", "y"));
            make_observable("TestGaussianResolution2D::NormalizationPDF(x,y)",
                            Unit::None(),
                            &TestGaussianResolution2DPDF::norm,
                            std::make_tuple("x_min", "x_max", "y_min", "y_max"));

            auto signal_pdf = make_signal_pdf("TestGaussianResolution2D::P(x,y)",
                                              "PDF for testing purpose only: 2D Gaussian resolution as a function of the offsets x, y.",
                                              Options{},
                                              "TestGaussianResolution2D::UnnormalizedPDF(x,y)",
                                              std::make_tuple("x"s, "y"s),
                                              "TestGaussianResolution2D::NormalizationPDF(x,y)",
                                              std::make_tuple("x_min"s, "x_max"s, "y_min"s, "y_max"s));
            SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
        }
    }
} // namespace eos::test
