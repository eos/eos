plot:
    title: 'Example Plot'
    x:
        label: '$q^2$'
        unit: '$\textnormal{GeV}^2$'
        range: [0.0, 26]
    y:
        label: '$\mathcal{B}(B\to \pi\ell\nu)$'
        unit: '$10^{-7}\,/\,\text{GeV}^2$'
        range: [1.0e-6, 1.0e-5]

contents:
    - name: 'B->pi mu nu (uncertainty)'
      type: 'uncertainty'
      color: 'red'
      opacity: 0.5
      hdf5-file: '/tmp/unc_btopilnu.hdf5'

    - name: 'B->pi mu nu'
      type: 'observable'
      color: 'red'
      observable: 'B->pilnu::dBR/ds;l=mu,model=SM'
      kinematic: 's'
      range: [0.01, 26.0]
      samples: 200

    - name: 'B->pi mu nu (V_ub inclusive)'
      type: 'observable'
      color: 'black'
      observable: 'B->pilnu::dBR/ds;l=mu,model=CKMScan'
      kinematic: 's'
      range: [0.01, 26.0]
      samples: 200
      parameters:
        'CKM::abs(V_ub)': 4.2e-3
