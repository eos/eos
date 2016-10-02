eos-propagate-uncertainty \
    --global-option model CKMScan \
    --global-option form-factors BCL2008 \
    --vary "CKM::abs(V_ub)"          2e-3   5e-3  --prior gaussian 4.19e-3 4.45e-3 4.71e-3 \
    --vary "decay-constant::B_u"     0.167  0.209 --prior gaussian 0.181   0.188   0.195   \
    --observable "B_u->lnu::BR;l=tau" \
    --workers 4 \
    --samples 100000 \
    --output /tmp/unc_btotaunu.hdf5
