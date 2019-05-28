plot:
    title: 'Example Plot'
    x:
        label: '$q^2$'
        unit: '$\textnormal{GeV}^2$'
        range: [0.0, 26]
    y:
        label: '$\mathcal{B}(B\to \pi\ell\nu)$'
        unit: '$10^{-5}\,/\,\textnormal{GeV}^2$'
        range: [1.0e-6, 1.0e-5]
        scale: 1e-5

contents:
    - name: 'B->pi mu nu (uncertainty)'
      type: 'uncertainty'
      color: 'red'
      opacity: 0.5
      hdf5-file: '/tmp/unc_btopilnu.hdf5'

    - name: 'B->pi mu nu'
      type: 'observable'
      color: 'red'
      observable: 'B->pilnu::dBR/dq2;l=mu,model=SM'
      kinematic: 'q2'
      range: [0.01, 26.0]
      samples: 200
      parameters-from-file: '/tmp/mode-btopipi-ff.yaml'

    - name: 'B->pi mu nu (V_ub inclusive)'
      type: 'observable'
      color: 'black'
      observable: 'B->pilnu::dBR/dq2;l=mu,model=CKMScan'
      kinematic: 'q2'
      range: [0.01, 26.0]
      samples: 200
      parameters:
        'CKM::abs(V_ub)':         4.20e-3
        'B->pi::f_+(0)@BCL2008':  2.66e-1
        'B->pi::b_+^1@BCL2008':  -2.67e+0
        'B->pi::b_+^2@BCL2008':  +2.23e-2

    - name: 'B->pi l nu (BaBar measurement)'
      type: 'constraint'
      constraints:
        - 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar-2010A'
        - 'B^0->pi^+lnu::BR[4.0,8.0]@BaBar-2010A'
        - 'B^0->pi^+lnu::BR[8.0,12.0]@BaBar-2010A'
      variable: 'q2'
      color: 'gray'
      rescale-by-width: true

    - type: 'watermark'
