/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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

#include <test/test.hh>
#include <eos/models/model.hh>
#include <eos/models/ckm.hh>
#include <eos/models/standard-model.hh>

#include <cmath>

using namespace test;
using namespace eos;

Parameters
reference_parameters()
{
    Parameters result = Parameters::Defaults();
    result["CKM::abs(V_ud)"] = 0.12345;
    result["CKM::abs(V_us)"] = 0.23456;
    result["CKM::abs(V_ub)"] = 0.34567;
    result["CKM::abs(V_cd)"] = 0.45678;
    result["CKM::abs(V_cs)"] = 0.56789;
    result["CKM::abs(V_cb)"] = 0.67890;
    result["CKM::abs(V_td)"] = 0.78901;
    result["CKM::abs(V_ts)"] = 0.89012;
    result["CKM::abs(V_tb)"] = 0.90123;
    result["CKM::arg(V_ud)"] = 1.12345;
    result["CKM::arg(V_us)"] = 1.23456;
    result["CKM::arg(V_ub)"] = 1.34567;
    result["CKM::arg(V_cd)"] = 1.45678;
    result["CKM::arg(V_cs)"] = 1.56789;
    result["CKM::arg(V_cb)"] = 1.67890;
    result["CKM::arg(V_td)"] = 1.78901;
    result["CKM::arg(V_ts)"] = 1.89012;
    result["CKM::arg(V_tb)"] = 1.90123;

    return result;
}

class CKMScanModelMakeTest :
    public TestCase
{
    public:
        CKMScanModelMakeTest() :
            TestCase("ckm_scan_model_make_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::shared_ptr<Model> m = Model::make("CKM", reference_parameters(), Options());
            }
            catch (NoSuchModelError &)
            {
                TEST_CHECK_FAILED("Model::make does not know the model 'CKM'");
            }
            catch (...)
            {
                throw;
                TEST_CHECK_FAILED("Unknown Exception while making 'CKM'");
            }
        }
} ckm_scan_model_make_test;

class CKMMatrixElementsTest :
    public TestCase
{
    public:
        CKMMatrixElementsTest() :
            TestCase("ckm_matrix_elements_test")
        {
        }

        virtual void run() const
        {
            /* Test passing of CKM matrix elements via polar parametrisations */
            {
                static const double eps = 1e-8;

                Parameters p = reference_parameters();
                Options o{ };

                CKMScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ud()), 0.12345, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_us()), 0.23456, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ub()), 0.34567, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cd()), 0.45678, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cs()), 0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cb()), 0.67890, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_td()), 0.78901, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ts()), 0.89012, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_tb()), 0.90123, eps);

                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_ud()), 1.12345, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_us()), 1.23456, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_ub()), 1.34567, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_cd()), 1.45678, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_cs()), 1.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_cb()), 1.67890, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_td()), 1.78901, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_ts()), 1.89012, eps);
                TEST_CHECK_NEARLY_EQUAL(arg(model.ckm_tb()), 1.90123, eps);
            }
        }
} ckm_matrix_elements_test;

class CKMDefaultValuesTest :
    public TestCase
{
    public:
        CKMDefaultValuesTest() :
            TestCase("ckm_default_values_test")
        {
        }

        virtual void run() const
        {
            /* Test passing of CKM matrix elements via polar parametrisations */
            {
                static const double eps = 1e-8;

                Parameters p = Parameters::Defaults();
                Options o{ };

                StandardModel sm(p);
                CKMScanModel ckm(p, o);

                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_ud()),    abs(ckm.ckm_ud()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_us()),    abs(ckm.ckm_us()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_ub()),    abs(ckm.ckm_ub()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_cd()),    abs(ckm.ckm_cd()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_cs()),    abs(ckm.ckm_cs()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_cb()),    abs(ckm.ckm_cb()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_td()),    abs(ckm.ckm_td()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_ts()),    abs(ckm.ckm_ts()),  eps);
                TEST_CHECK_NEARLY_EQUAL( abs(sm.ckm_tb()),    abs(ckm.ckm_tb()),  eps);

                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_ud()),    arg(ckm.ckm_ud()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_us()),    arg(ckm.ckm_us()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_ub()),    arg(ckm.ckm_ub()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_cd()),    arg(ckm.ckm_cd()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_cs()),    arg(ckm.ckm_cs()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_cb()),    arg(ckm.ckm_cb()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_td()),    arg(ckm.ckm_td()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_ts()),    arg(ckm.ckm_ts()),  eps);
                TEST_CHECK_NEARLY_EQUAL( arg(sm.ckm_tb()),    arg(ckm.ckm_tb()),  eps);
            }
        }
} ckm_default_values_test;
