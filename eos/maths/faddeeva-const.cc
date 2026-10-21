/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Tables for Faddeeva's scaled complex error function, see faddeeva.hh.
 *
 * Translated to C++ from the library libcerf,
 *   https://jugit.fz-juelich.de/mlz/lib/cerf
 *
 * Copyright:
 *   (C) 2012      Massachusetts Institute of Technology
 *   (C) 2013-2026 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, 2025
 */

#include <eos/maths/faddeeva-const.hh>

namespace eos
{
    namespace faddeeva_impl
    {
        namespace imwofx
        {
            // clang-format off
            alignas(64) const double coeffs0[144 * 8] = {
                -0x1.5f64e22c62f72p-77, -0x1.4966ce7e468c4p-70, 0x1.05741c6266c93p-61, 0x1.dba728f22ea3ap-55, -0x1.23a3c180da8aep-45, -0x1.dfffcd55d84a5p-41, 0x1.c4b1a0a360ce9p-30, -0x1.35c476fb1ab4ap-24, // subdomain 0:0 (0.5..0.515625)
                -0x1.4486166055455p-77, -0x1.63d0423100db9p-70, 0x1.f2cd9ddc0065ap-62, 0x1.0dc38852202ddp-54, -0x1.1cb24fa01b080p-45, -0x1.5c170679a3826p-40, 0x1.c1d1f91a89952p-30, -0x1.fab2d97783a0fp-25, // subdomain 0:1 (0.515625..0.53125)
                -0x1.285c8940e85e8p-77, -0x1.7c1315a102464p-70, 0x1.d8eac748f5a48p-62, 0x1.2c23614adec51p-54, -0x1.14e6b5d0236e5p-45, -0x1.c56a5c8f207b6p-40, 0x1.bde77b72f3e02p-30, -0x1.8ab62e7cf5347p-25, // subdomain 0:2 (0.53125..0.546875)
                -0x1.0b178b3faff13p-77, -0x1.921748ee304a8p-70, 0x1.bd674ac5411b8p-62, 0x1.48d7d6466db52p-54, -0x1.0c4c45be76c70p-45, -0x1.15d519a50fc8ap-39, 0x1.b8f975aa5faa2p-30, -0x1.1bd4bcf55cbf1p-25, // subdomain 0:3 (0.546875..0.5625)
                -0x1.d9ced5d3b9de6p-78, -0x1.a5c898ff0d1b7p-70, 0x1.a06c4c87ef58ap-62, 0x1.63c82449aaacap-54, -0x1.02ef075173968p-45, -0x1.47459b313965dp-39, 0x1.b30ff7992d66ep-30, -0x1.5c9cd1b4c00e6p-26, // subdomain 0:4 (0.5625..0.578125)
                -0x1.9bfa39b274eb1p-78, -0x1.b7168fea57eccp-70, 0x1.82243b5fd1a56p-62, 0x1.7cde2545496bcp-54, -0x1.f1b74d565506bp-46, -0x1.76e359f1ec457p-39, 0x1.ac33c7ac5def9p-30, -0x1.0983ec7cb5bb3p-27, // subdomain 0:5 (0.578125..0.59375)
                -0x1.5d13c7e82b73ap-78, -0x1.c5f48e8177409p-70, 0x1.62ba8be44068ep-62, 0x1.94066299a7fe6p-54, -0x1.dc3ec33737a02p-46, -0x1.a48d6b26eb5e6p-39, 0x1.a46e5705a5c50p-30, 0x1.3dc053ddacb64p-28, // subdomain 0:6 (0.59375..0.609375)
                -0x1.1d7df550fdb74p-78, -0x1.d259cf089b1e4p-70, 0x1.425b7337a8ee6p-62, 0x1.a9302349a75c4p-54, -0x1.c58fe940e2e72p-46, -0x1.d0256bd4015e5p-39, 0x1.9bc9b512f6335p-30, 0x1.1f872d166753fp-26, // subdomain 0:7 (0.609375..0.625)
                -0x1.bb34f728d0759p-79, -0x1.dc416136ccddap-70, 0x1.2133a1d5da961p-62, 0x1.bc4d75d97705cp-54, -0x1.adc70f20f7c10p-46, -0x1.f98f96b2b88afp-39, 0x1.925082af4a3b1p-30, 0x1.eb165be633051p-26, // subdomain 0:8 (0.625..0.640625)
                -0x1.3b9366d03dab7p-79, -0x1.e3aa1fa4f5ff9p-70, 0x1.fedffdcb76be1p-63, 0x1.cd5335ddd02ecp-54, -0x1.9501342264fffp-46, -0x1.10596b3f58d4bp-38, 0x1.880de4ded41e5p-30, 0x1.58db0e9ac6054p-25, // subdomain 0:9 (0.640625..0.65625)
                -0x1.79a7d230660cep-80, -0x1.e8969edf2ce45p-70, 0x1.ba7ac905bb18ap-63, 0x1.dc390d477f24ap-54, -0x1.7b5bdfe87f0bbp-46, -0x1.22bc6a4aa6516p-38, 0x1.7d0d7738db0f9p-30, 0x1.b9824ec639c1cp-25, // subdomain 0:10 (0.65625..0.671875)
                -0x1.fd6fca01ebc7cp-82, -0x1.eb0d1657346b1p-70, 0x1.7590baded5ac9p-63, 0x1.e8f9717e420d5p-54, -0x1.60f4fb3564f87p-46, -0x1.33e701f376d5ep-38, 0x1.715b3e11aeff3p-30, 0x1.0ba97c809f39ap-24, // subdomain 0:11 (0.671875..0.6875)
                0x1.d96a220f36b7ep-82, -0x1.eb174381fbf51p-70, 0x1.3079c93413559p-63, 0x1.f3919c64e5d50p-54, -0x1.45eaa8f5aedd2p-46, -0x1.43d0d41ed46a8p-38, 0x1.650398771344cp-30, 0x1.391112477c444p-24, // subdomain 0:12 (0.6875..0.703125)
                0x1.6628252760084p-80, -0x1.e8c24773fb718p-70, 0x1.d7189b7752352p-64, 0x1.fc01816718e08p-54, -0x1.2a5b1fc9c82e2p-46, -0x1.5272ecf13d889p-38, 0x1.5813321160758p-30, 0x1.64e4052a64b3ap-24, // subdomain 0:13 (0.703125..0.71875)
                0x1.276519d1e3c77p-79, -0x1.e41e7f5970d14p-70, 0x1.4e391d7108830p-64, 0x1.0125df5bcdcd7p-53, -0x1.0e64844366921p-46, -0x1.5fc7c2d7b5282p-38, 0x1.4a96f4fb62e58p-30, 0x1.8f100b1df933cp-24, // subdomain 0:14 (0.71875..0.734375)
                0x1.97884f6c72882p-79, -0x1.dd3f582fb9dd9p-70, 0x1.8df13650335abp-65, 0x1.033ac5759f704p-53, -0x1.e44988168ef79p-47, -0x1.6bcb34d0f31a8p-38, 0x1.3c9bfba2a6fa5p-30, 0x1.b7847a84efeb6p-24, // subdomain 0:15 (0.734375..0.75)
                0x1.0177e415d445ap-78, -0x1.d43b1e2c39c7bp-70, 0x1.07d2e00951b00p-66, 0x1.04434f913b35dp-53, -0x1.ab72e44195a1fp-47, -0x1.767a8702c4416p-38, 0x1.2e2f82c177ad6p-30, 0x1.de324f59b144ep-24, // subdomain 0:16 (0.75..0.765625)
                0x1.348c395f5c7f3p-78, -0x1.c92ac84362650p-70, -0x1.00ee3e1ec8a76p-66, 0x1.04448efc10d15p-53, -0x1.727f48c44aed9p-47, -0x1.7fd45da9ec5cap-38, 0x1.1f5edb834f9b5p-30, 0x1.0186174cd7e39p-23, // subdomain 0:17 (0.765625..0.78125)
                0x1.64c55cf5ae53ap-78, -0x1.bc29c0478e1e8p-70, -0x1.7e113c55adae1p-65, 0x1.0344bbb844470p-53, -0x1.39a7a51195449p-47, -0x1.87d8b674f9820p-38, 0x1.10375de3e3132p-30, 0x1.130333fc66b72p-23, // subdomain 0:18 (0.78125..0.796875)
                0x1.91ed86d75a529p-78, -0x1.ad55a81b5400ep-70, -0x1.39eadde59eef1p-64, 0x1.014b2513afcd1p-53, -0x1.01236d84aaf09p-47, -0x1.8e88e05c70f76p-38, 0x1.00c65b583b41bp-30, 0x1.238b7af28d87fp-23, // subdomain 0:19 (0.796875..0.8125)
                0x1.bbd5296c28196p-78, -0x1.9cce1c75e0b34p-70, -0x1.b0635737a3233p-64, 0x1.fcc041feeea32p-54, -0x1.9250c02ad5667p-48, -0x1.93e7720c8ba07p-38, 0x1.e23223a132b79p-31, 0x1.331abdf290bc9p-23, // subdomain 0:20 (0.8125..0.828125)
                0x1.e2531c704a355p-78, -0x1.8ab475b9a58fep-70, -0x1.10fe1b6ba150dp-63, 0x1.f519f49533553p-54, -0x1.23d499a2b7886p-48, -0x1.97f83ef65b48dp-38, 0x1.c2793e420c6fep-31, 0x1.41ada1ae56e95p-23, // subdomain 0:21 (0.828125..0.84375)
                0x1.02a25d7cd5235p-77, -0x1.772b875e5d8cbp-70, -0x1.4724a0b674ddbp-63, 0x1.ebb7bbc602456p-54, -0x1.6e6b894f3fb51p-49, -0x1.9ac04b1fa80a7p-38, 0x1.a27be9bb21acfp-31, 0x1.4f419b62e4b75p-23, // subdomain 0:22 (0.84375..0.859375)
                0x1.1246fa576ac01p-77, -0x1.62575e6f034fep-70, -0x1.7a741a7fbc3bfp-63, 0x1.e0af8e82c4664p-54, -0x1.334bf37e32f85p-50, -0x1.9c45bdca1c713p-38, 0x1.8253996c9d676p-31, 0x1.5bd4edb68b4e5p-23, // subdomain 0:23 (0.859375..0.875)
                0x1.200ca930a81aap-77, -0x1.4c5cff9ae471dp-70, -0x1.aac08bf0aa5aep-63, 0x1.d418d8497831dp-54, 0x1.af9f993a67f5fp-52, -0x1.9c8fd30b6ebf4p-38, 0x1.62193abc6a7b7p-31, 0x1.6766a4e11f6efp-23, // subdomain 0:24 (0.875..0.890625)
                0x1.2bebf72b1ddbbp-77, -0x1.356225568d773p-70, -0x1.d7e33a8893719p-63, 0x1.c60c4f8518954p-54, 0x1.ff727c27d292ap-50, -0x1.9ba6cc711d2c2p-38, 0x1.41e5200dbc92cp-31, 0x1.71f6922738c16p-23, // subdomain 0:25 (0.890625..0.90625)
                0x1.35e0b050085abp-77, -0x1.1d8cfe85e4841p-70, -0x1.00dd6188f8548p-62, 0x1.b6a3cb217a20fp-54, 0x1.c30bf83080a55p-49, -0x1.9993e0cb09cd3p-38, 0x1.21ceecf4779c4p-31, 0x1.7b8546b211961p-23, // subdomain 0:26 (0.90625..0.921875)
                0x1.3de9d001454a0p-77, -0x1.0503ee247c22dp-70, -0x1.141592bbcf6ccp-62, 0x1.a5fa179c150c2p-54, 0x1.3facdbfc136b4p-48, -0x1.96612b38c9196p-38, 0x1.01ed83b456d0dp-31, 0x1.84140dce27046p-23, // subdomain 0:27 (0.921875..0.9375)
                0x1.44096bf749d53p-77, -0x1.d7da98b806e81p-71, -0x1.258ee3d304de0p-62, 0x1.942acbda150b0p-54, 0x1.9a0dfd9b50037p-48, -0x1.92199995cc0eap-38, 0x1.c4ade833083f2p-32, 0x1.8ba4e69b2f9eep-23, // subdomain 0:28 (0.9375..0.953125)
                0x1.48449a1b75028p-77, -0x1.a4de52e59dd57p-71, -0x1.3540b7d9a4a61p-62, 0x1.81521e0a6c9a3p-54, 0x1.f06c3b12bf7ffp-48, -0x1.8cc8da60b5a00p-38, 0x1.8640d76c37a9dp-32, 0x1.923a7d3863a81p-23, // subdomain 0:29 (0.953125..0.96875)
                0x1.4aa351979863cp-77, -0x1.715e25ebae037p-71, -0x1.43251acbd4d31p-62, 0x1.6d8cb8d896d5fp-54, 0x1.214853eed7062p-47, -0x1.867b4a3a234e3p-38, 0x1.48bc4f1bc4334p-32, 0x1.97d8237864ba6p-23, // subdomain 0:30 (0.96875..0.984375)
                0x1.4b30477b0ccccp-77, -0x1.3da3b6f4998afp-71, -0x1.4f38b666c7425p-62, 0x1.58f7913230b0cp-54, 0x1.482552d7c0df9p-47, -0x1.7f3de116fbfc0p-38, 0x1.0c46ce6e781d5p-32, 0x1.9c81c928542ecp-23, // subdomain 0:31 (0.984375..1.0)
                0x1.48b0e0d084a75p-67, -0x1.e07b368c13c0ep-63, -0x1.5df744fb97295p-54, 0x1.38d1e7cfc8adep-47, 0x1.7e1fad6d94d1fp-41, -0x1.72bdb5c1a62e6p-33, 0x1.67c3ab27918b2p-29, 0x1.a1c0cbf7896d5p-20, // subdomain 1:0 (1.0..1.03125)
                0x1.3f987df64222bp-67, -0x1.15903ec2d6504p-63, -0x1.6b457fb24121ep-54, 0x1.0c2ca86dcafb3p-47, 0x1.bdb347dc7be9bp-41, -0x1.5f4cff0b01180p-33, 0x1.0bf12add8f6f0p-30, 0x1.a598373fc1d60p-20, // subdomain 1:1 (1.03125..1.0625)
                0x1.3071aac8c9976p-67, -0x1.4907ff426ef03p-65, -0x1.719069f2fa83cp-54, 0x1.bc9cfd06cb0c9p-48, 0x1.f35c6b4e01a0ep-41, -0x1.491c9775a5ae9p-33, -0x1.3b27cf13c5262p-31, 0x1.a601f4ad19298p-20, // subdomain 1:2 (1.0625..1.09375)
                0x1.1bf7b9f5fe62ep-67, 0x1.977f19f433480p-65, -0x1.712cf37830b5cp-54, 0x1.6025a712a4550p-48, 0x1.0f776ac7de11bp-40, -0x1.30a4b62560955p-33, -0x1.14f197806a9e5p-29, 0x1.a335643107027p-20, // subdomain 1:3 (1.09375..1.125)
                0x1.02fb2b795e61ap-67, 0x1.0fc1604f77656p-63, -0x1.6a8df715ec92cp-54, 0x1.0491134681023p-48, 0x1.20367b14723bap-40, -0x1.165e8a6fe71e9p-33, -0x1.cb378a07b0abap-29, 0x1.9d6f9839f90d5p-20, // subdomain 1:4 (1.125..1.15625)
                0x1.ccb590803b4a0p-68, 0x1.a8d63b5cd9ad1p-63, -0x1.5e3fa8a8f79eep-54, 0x1.56ba19a73e94cp-49, 0x1.2c0266a186d3dp-40, -0x1.f5843b1a00478p-34, -0x1.38503e43a2ae1p-28, 0x1.94f2272e1047ep-20, // subdomain 1:5 (1.15625..1.1875)
                0x1.8df9997b04ad6p-68, 0x1.178fd7d500816p-62, -0x1.4ce2b2a5a7238p-54, 0x1.5786aadadd5d0p-50, 0x1.33060cab33af5p-40, -0x1.bc88bc38bfdeep-34, -0x1.823da33bb81f3p-28, 0x1.8a01ffbc6dbe1p-20, // subdomain 1:6 (1.1875..1.21875)
                0x1.4b90ac08ef682p-68, 0x1.5097ca970e478p-62, -0x1.372729a23dedcp-54, 0x1.5334cbd98b4f3p-54, 0x1.357dfad92f421p-40, -0x1.82aaba0a47150p-34, -0x1.c32f32a2e047ep-28, 0x1.7ce644e6af177p-20, // subdomain 1:7 (1.21875..1.25)
                0x1.073e7ea063456p-68, 0x1.7f136f4bedb51p-62, -0x1.1dc77ac4eb0dcp-54, -0x1.1583f160ba01dp-50, 0x1.33b606913ca2ep-40, -0x1.48bd832376fc0p-34, -0x1.fb123984d8396p-28, 0x1.6de736293241ap-20, // subdomain 1:8 (1.25..1.28125)
                0x1.85625c4af9068p-69, 0x1.a2d8460838ea9p-62, -0x1.01837787b3ff6p-54, -0x1.12ad4c11e24a9p-49, 0x1.2e06c20a50d61p-40, -0x1.0f84f2047abcfp-34, -0x1.14f9ee0c9b817p-27, 0x1.5d4d279ef90ebp-20, // subdomain 1:9 (1.28125..1.3125)
                0x1.fdde57c4f0217p-70, 0x1.bbfcaf75da630p-62, -0x1.c6373be6a97fbp-55, -0x1.8be6589abcb52p-49, 0x1.24d2de6d692cap-40, -0x1.af6690a8ca2a6p-35, -0x1.27ff3f8f01066p-27, 0x1.4b5f8d69b5695p-20, // subdomain 1:10 (1.3125..1.34375)
                0x1.f7ca2290cc9c5p-71, 0x1.cad26428fc940p-62, -0x1.869970e07c853p-55, -0x1.f58a9e0f288d9p-49, 0x1.18848b7eb0f57p-40, -0x1.43cf10634c7dfp-35, -0x1.36bb72a5d5ad8p-27, 0x1.38641d0bbebe5p-20, // subdomain 1:11 (1.34375..1.375)
                0x1.3814b739e3f01p-75, 0x1.cfdfcc71a47c1p-62, -0x1.4597d7902fd0cp-55, -0x1.278a1e1ed3104p-48, 0x1.098ae418757fbp-40, -0x1.bab12a65e3c5dp-36, -0x1.415c8a45ba583p-27, 0x1.249e06c466d1fp-20, // subdomain 1:12 (1.375..1.40625)
                -0x1.a8df8f7b521fep-71, 0x1.cbd87f0ca7d0bp-62, -0x1.04868a3271ae3p-55, -0x1.4c2a9f40f6ec4p-48, 0x1.f0aeea828cc58p-41, -0x1.f3b8796c2cc93p-37, -0x1.481a8b1cfecdfp-27, 0x1.104d4863105f9p-20, // subdomain 1:13 (1.40625..1.4375)
                -0x1.9b60be9d5d74bp-70, 0x1.bf952df96a4a0p-62, -0x1.893b51f997a96p-56, -0x1.68b8a8853a00ep-48, 0x1.cab7d7d0cbd14p-41, -0x1.1a7ff32391089p-38, -0x1.4b35d27ca710ep-27, 0x1.f75c34dbe9acap-21, // subdomain 1:14 (1.4375..1.46875)
                -0x1.244f53d20532fp-69, 0x1.ac0b4747a9e61p-62, -0x1.0de970cd7ef9dp-56, -0x1.7d6b142e67676p-48, 0x1.a20fe6d193249p-41, 0x1.775e47454e5fep-38, -0x1.4af56cba81d81p-27, 0x1.cdf0efd36e51ap-21, // subdomain 1:15 (1.46875..1.5)
                -0x1.6d45cfe3acfbcp-69, 0x1.92448c26971edp-62, -0x1.31fcb99305c6fp-57, -0x1.8a99732205bc8p-48, 0x1.778e99d974f02p-41, 0x1.e561a1b68d52ep-37, -0x1.47a573075ca04p-27, 0x1.a4bf8112a2228p-21, // subdomain 1:16 (1.5..1.53125)
                -0x1.a82bc2456488ep-69, 0x1.7356db276ee4fp-62, -0x1.60a2b956bd09cp-59, -0x1.90b763d7b220ap-48, 0x1.4c004b76d24b4p-41, 0x1.7760f40edf348p-36, -0x1.419577dfbce2ap-27, 0x1.7c24e9297261bp-21, // subdomain 1:17 (1.53125..1.5625)
                -0x1.d4f97bea9de61p-69, 0x1.505c63d6b4c10p-62, 0x1.bc4aabacae315p-59, -0x1.904fb29fd6a72p-48, 0x1.202360dc52047p-41, 0x1.eba747629d6c8p-36, -0x1.391708523c8ffp-27, 0x1.547412d2f671ap-21, // subdomain 1:18 (1.5625..1.59375)
                -0x1.f3ff18a149bd8p-69, 0x1.2a6c71e2ef0a1p-62, 0x1.21bd0792680aap-57, -0x1.89ff66499d7e9p-48, 0x1.e94c0c4ab8335p-42, 0x1.27c3c0ff51961p-35, -0x1.2e7c475c6102bp-27, 0x1.2df5aa39194d4p-21, // subdomain 1:19 (1.59375..1.625)
                -0x1.02ed206105ab5p-68, 0x1.0294f55b6a810p-62, 0x1.be78ea8a3f9b4p-57, -0x1.7e70d989038f7p-48, 0x1.9448e99bb8ab3p-42, 0x1.51a1ace01fcdfp-35, -0x1.2216a79cdebabp-27, 0x1.08e81ef149cb7p-21, // subdomain 1:20 (1.625..1.65625)
                -0x1.05b55c8098d94p-68, 0x1.b3a9c90111c50p-63, 0x1.223e5fde88c8cp-56, -0x1.6e56fe2e21f96p-48, 0x1.424f82b1f012ep-42, 0x1.73a9c553eb3eap-35, -0x1.1435c697a2ca9p-27, 0x1.caff996222f2cp-22, // subdomain 1:21 (1.65625..1.6875)
                -0x1.02e31a0f95061p-68, 0x1.622f15004f265p-63, 0x1.59c48833126c4p-56, -0x1.5a68e36edc5a2p-48, 0x1.e894a8a3a05f8p-43, 0x1.8e3041eca0e34p-35, -0x1.05266bcaf9ef8p-27, 0x1.87ce8cc39b07fp-22, // subdomain 1:22 (1.6875..1.71875)
                -0x1.f62bd52217accp-69, 0x1.12619b35c22d3p-63, 0x1.85ee2c759972cp-56, -0x1.435d945d0a33dp-48, 0x1.55f77c023a48ap-43, 0x1.a19d748d17850p-35, -0x1.ea6359c6ad628p-28, 0x1.487f7e0438362p-22, // subdomain 1:23 (1.71875..1.75)
                -0x1.ddf7c87bbeaaep-69, 0x1.8b69bb5e98d67p-64, 0x1.a712863a1188dp-56, -0x1.29e85f4cba3c4p-48, 0x1.9bf9b4a7ed38ep-44, 0x1.ae6a0a0ffa310p-35, -0x1.c9386ef299269p-28, 0x1.0d4317505ea27p-22, // subdomain 1:24 (1.75..1.78125)
                -0x1.be974fa0e8124p-69, 0x1.f591eea2d13abp-65, 0x1.bdb928be93c29p-56, -0x1.0eb5928dcba1bp-48, 0x1.46345fd17f0d4p-45, 0x1.b51b5e5a2aad3p-35, -0x1.a74b7fee706e3p-28, 0x1.ac72e7df7b605p-23, // subdomain 1:25 (1.78125..1.8125)
                -0x1.997d8dc526ab7p-69, 0x1.d3934993e98d7p-66, 0x1.ca91fed11576fp-56, -0x1.e4cf6d27c4ca0p-49, -0x1.eb8cd44609571p-47, 0x1.b63ff4f642e70p-35, -0x1.8511291971411p-28, 0x1.46e6d9dd60e23p-23, // subdomain 1:26 (1.8125..1.84375)
                -0x1.7019564cde0a7p-69, -0x1.2c9645532d5ffp-70, 0x1.ce6d273874e0fp-56, -0x1.ab2a932b3af95p-49, -0x1.04f4c7786a7e9p-44, 0x1.b26c24c6216fap-35, -0x1.62f2a4ef5cf52p-28, 0x1.d3cff2d5a044fp-24, // subdomain 1:27 (1.84375..1.875)
                -0x1.43cc2eb897121p-69, -0x1.c373d6356026bp-66, 0x1.ca32e8e3cef4ep-56, -0x1.718df9d23e3a2p-49, -0x1.b33834ad978abp-44, 0x1.aa3702b160affp-35, -0x1.414d76cc53f88p-28, 0x1.2ac6c4f25b8d4p-24, // subdomain 1:28 (1.875..1.90625)
                -0x1.15e2a6e6baca2p-69, -0x1.9dcd9d1071579p-65, 0x1.bedbf167779eep-56, -0x1.38ece20904bacp-49, -0x1.243aeb6cfdd09p-43, 0x1.9e37948ff0a38p-35, -0x1.20735095efe93p-28, 0x1.24c12961e9c77p-25, // subdomain 1:29 (1.90625..1.9375)
                -0x1.cf1c3826e5674p-70, -0x1.1e7f80ddf7e00p-64, 0x1.ad6a08cf0a80ep-56, -0x1.021ad4374b07dp-49, -0x1.62a594ab56873p-43, 0x1.8f0252f87648fp-35, -0x1.00aa2d5e20d16p-28, 0x1.44bc72a28f3efp-29, // subdomain 1:30 (1.9375..1.96875)
                -0x1.73bfdadd75b19p-70, -0x1.5fb131a3f9846p-64, 0x1.96e15dd903bdep-56, -0x1.9b95d38dd129ep-50, -0x1.95558ed6c93b4p-43, 0x1.7d26fe21eb9f9p-35, -0x1.c45935505612ap-29, -0x1.ba03ac4b053dcp-26, // subdomain 1:31 (1.96875..2.0)
                -0x1.e2736d427bc21p-61, -0x1.a70e76e9b8c80p-55, 0x1.6dc6d15996459p-48, -0x1.0a7ecf5a68346p-43, -0x1.ccbf2d5c0f3a1p-37, 0x1.5e908bf4384f9p-30, -0x1.6e863bce988f8p-25, -0x1.07d168cd228fep-21, // subdomain 2:0 (2.0..2.0625)
                -0x1.5a2d1c64b07a6p-62, -0x1.d99a019cc69fbp-55, 0x1.2e11e5457258fp-48, -0x1.8d0e01644bc53p-45, -0x1.f42b9dc145e26p-37, 0x1.313c5a0461632p-30, -0x1.07f5a6b6158f4p-25, -0x1.a4d917faf2c60p-21, // subdomain 2:1 (2.0625..2.125)
                0x1.383e076449703p-63, -0x1.e036da68108e8p-55, 0x1.d4effc155e81dp-49, 0x1.173cec1d35847p-46, -0x1.fa97cbed30569p-37, 0x1.01d07d6480600p-30, -0x1.5fef6a5d748e1p-26, -0x1.091b1ace5aeefp-20, // subdomain 2:2 (2.125..2.1875)
                0x1.08146bc4e5d50p-61, -0x1.c4a2b6efa0cbdp-55, 0x1.50f63456088aap-49, 0x1.0efbffe36924bp-44, -0x1.e75a9caae42e5p-37, 0x1.a61a449529f74p-31, -0x1.9b04b6f2bd1d1p-27, -0x1.2ba4383fc198cp-20, // subdomain 2:3 (2.1875..2.25)
                0x1.816baa870e403p-61, -0x1.91172af8ff15dp-55, 0x1.b0c738fb92f07p-50, 0x1.98c3d182a795dp-44, -0x1.c1a163ba86363p-37, 0x1.4e066a7a4dcd1p-31, -0x1.5e322b5bea965p-28, -0x1.3dabb6bbb0161p-20, // subdomain 2:4 (2.25..2.3125)
                0x1.c01984de37486p-61, -0x1.4f4623fd0a770p-55, 0x1.c248ab1279bd8p-51, 0x1.ea4b7f63f44f6p-44, -0x1.8ffd65d31658fp-37, 0x1.fc77da796cd5ap-32, 0x1.07e9e55660532p-32, -0x1.42a000eefbd67p-20, // subdomain 2:5 (2.3125..2.375)
                0x1.cdc579a9c311fp-61, -0x1.07a1dc5f15eafp-55, 0x1.c4ac706423188p-53, 0x1.065e8d3827a5ap-43, -0x1.58189856b2350p-37, 0x1.70d5bc695506bp-32, 0x1.218f0322e6f5ap-28, -0x1.3d9d99b3710b6p-20, // subdomain 2:6 (2.375..2.4375)
                0x1.b5e914a9deb0fp-61, -0x1.81e0c24cd314bp-56, -0x1.1e7c678367ca7p-52, 0x1.050c872899743p-43, -0x1.1e8e1c0c2379ep-37, 0x1.f531dd36d4008p-33, 0x1.e20199627860fp-28, -0x1.315e1376a6eefp-20, // subdomain 2:7 (2.4375..2.5)
                0x1.8445d5d81e08cp-61, -0x1.00454b8fb5432p-56, -0x1.42ff625109fd6p-51, 0x1.ec42b64c465bfp-44, -0x1.cdc160dfa3e7cp-38, 0x1.335cca75de465p-33, 0x1.2fa0a60bea5a8p-27, -0x1.2030a858829cep-20, // subdomain 2:8 (2.5..2.5625)
                0x1.43b9fcc0d0324p-61, -0x1.2187cafd1449fp-57, -0x1.b2af682bef2c1p-51, 0x1.bc4e53a60f3fcp-44, -0x1.67106ec75a3b2p-38, 0x1.3393314dd862bp-34, 0x1.532738586c762p-27, -0x1.0bfa6e3da156ep-20, // subdomain 2:9 (2.5625..2.625)
                0x1.faf65156e647fp-62, -0x1.b441e3e30bff5p-59, -0x1.e9897a38c1c11p-51, 0x1.820ead3a19e79p-44, -0x1.0c1e56abfa716p-38, 0x1.24b78bcb6a45ap-36, 0x1.61950c8c2075fp-27, -0x1.ec784b963a570p-21, // subdomain 2:10 (2.625..2.6875)
                0x1.7160e744872e0p-62, 0x1.ba400620db51ap-61, -0x1.f3cfdde972e92p-51, 0x1.43db2b6eb0814p-44, -0x1.7ce5aeb595abbp-39, -0x1.87ba4edde3027p-36, 0x1.603c6dcca7eb2p-27, -0x1.c03797190f80fp-21, // subdomain 2:11 (2.6875..2.75)
                0x1.e972c78c5cdafp-63, 0x1.ecaaf29fdaa48p-59, -0x1.dd8a4b051df37p-51, 0x1.0683fa6f10723p-44, -0x1.f975f87c0ad51p-40, -0x1.afcba30a99493p-35, 0x1.53aa0a6e8353cp-27, -0x1.94e0a5082cb0ap-21, // subdomain 2:12 (2.75..2.8125)
                0x1.16131c1ae5ea7p-63, 0x1.6d302ccf377aep-58, -0x1.b1a3ee5c536c5p-51, 0x1.9ad0c302f92d6p-45, -0x1.2d2538a968aa0p-40, -0x1.22bc5f1027130p-34, 0x1.3f9fe1b8fce2ap-27, -0x1.6b9c773e1691dp-21, // subdomain 2:13 (2.8125..2.875)
                0x1.b7cf7ba40b5bdp-65, 0x1.a8af3f724af19p-58, -0x1.796294175764ep-51, 0x1.35475b5d57757p-45, -0x1.2046c5cc36243p-41, -0x1.4bc96dff751f5p-34, 0x1.271d9f5c0223dp-27, -0x1.45281931efe67p-21, // subdomain 2:14 (2.875..2.9375)
                -0x1.db12ba0acfcacp-68, 0x1.b676f4a7ea7b7p-58, -0x1.3c29a032da2f1p-51, 0x1.bd194fa861010p-46, -0x1.c5bfe51bc49e9p-44, -0x1.5b5d94f468deap-34, 0x1.0c71f3dfd4196p-27, -0x1.21ebe773e7231p-21, // subdomain 2:15 (2.9375..3.0)
                -0x1.8faf65dcfaafdp-65, 0x1.a3944f584e4e5p-58, -0x1.feeffaf4d2198p-52, 0x1.2e4b6679620dap-46, 0x1.a6aa806f3e1c0p-43, -0x1.589546155c18dp-34, 0x1.e2a24189fdc45p-28, -0x1.02104fc3e804ep-21, // subdomain 2:16 (3.0..3.0625)
                -0x1.2e2fc03e6b331p-64, 0x1.7b92a5a4bb8b1p-58, -0x1.8e1ee90c278a2p-52, 0x1.79c1f66387a78p-47, 0x1.a829f94c9f1fep-42, -0x1.494d180002450p-34, 0x1.addb6169c9393p-28, -0x1.cb1f2baf9d2b6p-22, // subdomain 2:17 (3.0625..3.125)
                -0x1.5e5eb7d153979p-64, 0x1.47fc7bb8d1cfep-58, -0x1.2a6a5a7071cb1p-52, 0x1.985c5b064637ep-48, 0x1.12d28e59064b2p-41, -0x1.322e5a2bc9989p-34, 0x1.7c222fcc5c889p-28, -0x1.9888fb1318e93p-22, // subdomain 2:18 (3.125..3.1875)
                -0x1.662ed369471b5p-64, 0x1.103277ad3fb73p-58, -0x1.abf6cdd58a996p-53, 0x1.32ba8c7adb965p-49, 0x1.30c5b05128625p-41, -0x1.16cb7dd708a23p-34, 0x1.4e655b41e321cp-28, -0x1.6bebf03795549p-22, // subdomain 2:19 (3.1875..3.25)
                -0x1.524766d9bbfecp-64, 0x1.b30228e368537p-59, -0x1.22681ab33489bp-53, -0x1.90d05197b0a48p-52, 0x1.37299d9d0e1ebp-41, -0x1.f38a9a157d10ap-35, 0x1.2518163551358p-28, -0x1.44c03448b1d4ep-22, // subdomain 2:20 (3.25..3.3125)
                -0x1.2d6b7ef0d9b80p-64, 0x1.4ebd8de10691ep-59, -0x1.6cf3d5e7ce2b3p-54, -0x1.1c62e1d459dadp-49, 0x1.2d8652135211bp-41, -0x1.b9e6f584e16e9p-35, 0x1.0054b17ea7ba2p-28, -0x1.22756c0629c33p-22, // subdomain 2:21 (3.3125..3.375)
                -0x1.00268023e5934p-64, 0x1.eefda38119e0fp-60, -0x1.93f209e5c468ep-55, -0x1.a83df928cf1a2p-49, 0x1.19c6e62bab304p-41, -0x1.831936c9aa1e8p-35, 0x1.bff1c30a9b910p-29, -0x1.047bfef7f44dep-22, // subdomain 2:22 (3.375..3.4375)
                -0x1.a195379be9f4ap-65, 0x1.5db7a13ed669dp-60, -0x1.4ffad051756f3p-56, -0x1.ee2c6cde798aep-49, 0x1.0065384c5b4b9p-41, -0x1.5095da76a85bfp-35, 0x1.8777fb561bea2p-29, -0x1.d4966f5f14dcep-23, // subdomain 2:23 (3.4375..3.5)
                -0x1.474b374525238p-65, 0x1.d30d022295f73p-61, -0x1.5faee86c2e779p-61, -0x1.014e50b5c2d90p-48, 0x1.c940c08b6d64cp-42, -0x1.2318ea8e90a2ep-35, 0x1.567fb4e30e2dcp-29, -0x1.a6c9e8c7b4945p-23, // subdomain 2:24 (3.5..3.5625)
                -0x1.ed50409ce4364p-66, 0x1.204656b6e2dd2p-61, 0x1.8d620d6e6380ap-57, -0x1.f5ebe9bdb088dp-49, 0x1.916c36abddacdp-42, -0x1.f5b74239cb3ddp-36, 0x1.2c3a103539c39p-29, -0x1.7eaf106d828ebp-23, // subdomain 2:25 (3.5625..3.625)
                -0x1.647fe2049a5c3p-66, 0x1.384d955603bd0p-62, 0x1.4159ddf9652a1p-56, -0x1.d4b5ed52f7f55p-49, 0x1.5c33a8f1d6166p-42, -0x1.af7eddd8f4f8ep-36, 0x1.07d696db786d7p-29, -0x1.5b7ca6755ded7p-23, // subdomain 2:26 (3.625..3.6875)
                -0x1.eaeb1a45e2e5ap-67, 0x1.f601bf3a3ff8ap-64, 0x1.7d05bf4ad188bp-56, -0x1.a84fe30909e0fp-49, 0x1.2b4bff2d595d9p-42, -0x1.72de5eb6736d0p-36, 0x1.d11c77effae85p-30, -0x1.3c82f9e80d14ap-23, // subdomain 2:27 (3.6875..3.75)
                -0x1.3de45ecea370dp-67, 0x1.4c18f23dfd27dp-72, 0x1.8d3da5603e8f2p-56, -0x1.7755f2376692ep-49, 0x1.ff187c955a9b8p-43, -0x1.3eef8e522eb81p-36, 0x1.9b53b9d90e15ap-30, -0x1.212a493579947p-23, // subdomain 2:28 (3.75..3.8125)
                -0x1.78d2279bd4492p-68, -0x1.37025c50069a9p-64, 0x1.813f920408aa0p-56, -0x1.4638ba3a3c6a5p-49, 0x1.b263f229f1a89p-43, -0x1.12b33ae8a571cp-36, 0x1.6d0afa1a28bd6p-30, -0x1.08f089ec6ab2ep-23, // subdomain 2:29 (3.8125..3.875)
                -0x1.7e320dbedfb6ep-69, -0x1.e5183099c6159p-64, 0x1.6477fa736eadfp-56, -0x1.17bf918719611p-49, 0x1.702da038f9bbfp-43, -0x1.da517fa9b04fdp-37, 0x1.45227a6f2fd40p-30, -0x1.e6cddf3c1c7c4p-24, // subdomain 2:30 (3.875..3.9375)
                -0x1.00c98e04fac4ep-70, -0x1.191b451bfa3d7p-63, 0x1.3f1bca8cc59e0p-56, -0x1.daf285f09351fp-50, 0x1.37b13caced3a1p-43, -0x1.9ab99aead5909p-37, 0x1.22a1b0e70124cp-30, -0x1.c05edaa100655p-24, // subdomain 2:31 (3.9375..4.0)
                0x1.64489b122c7eap-61, -0x1.1a22f4ed73d77p-54, 0x1.02abaed58347bp-48, -0x1.6e9a885608629p-43, 0x1.e5f0479d9ecc1p-38, -0x1.4d24704eca605p-32, 0x1.ee77ea69d2085p-27, -0x1.8e1666fab6a3ep-21, // subdomain 3:0 (4.0..4.125)
                0x1.7efb4c46c550ep-60, -0x1.d48af8a109457p-55, 0x1.7261b2334a4bfp-49, -0x1.002b776fa0d6ap-43, 0x1.5f0c48b74d092p-38, -0x1.fedba443b111bp-33, 0x1.9329572605ea7p-27, -0x1.563cf5917c5c6p-21, // subdomain 3:1 (4.125..4.25)
                0x1.7818fda642b05p-60, -0x1.5abb469d8b45dp-55, 0x1.ff9d9b49c75d7p-50, -0x1.651d5133ec62bp-44, 0x1.00fee147897e3p-38, -0x1.8de8b27011568p-33, 0x1.4c9f4e4e4cd10p-27, -0x1.286f648b42289p-21, // subdomain 3:2 (4.25..4.375)
                0x1.29252b97b611ep-60, -0x1.e1937d79bdc1bp-56, 0x1.5bbcbf835b2ebp-50, -0x1.f5ac60d000e2dp-45, 0x1.7e7e8731e65bfp-39, -0x1.3a969a9c2abe6p-33, 0x1.15475c12b48d3p-27, -0x1.0273a1fcaa6a5p-21, // subdomain 3:3 (4.375..4.5)
                0x1.a68a0da7bcfa5p-61, -0x1.43352523cb2ddp-56, 0x1.d811f1762d45bp-51, -0x1.656de48d77902p-45, 0x1.218bec98e2509p-39, -0x1.f8270a0b2444dp-34, 0x1.d278116f87a2cp-28, -0x1.c54ada51f47dcp-22, // subdomain 3:4 (4.5..4.625)
                0x1.1b8f5faa30310p-61, -0x1.abda57e412304p-57, 0x1.435fd170b1b5fp-51, -0x1.032b09e8f296cp-45, 0x1.bd90c28ac81b5p-40, -0x1.98c2c0630cdadp-34, 0x1.8b7916ca81d27p-28, -0x1.8f9378a62ce2ep-22, // subdomain 3:5 (4.625..4.75)
                0x1.7184d040c595ep-62, -0x1.1b734cda71b2ap-57, 0x1.c2103db091d0dp-52, -0x1.7ec70cd323f3dp-46, 0x1.5be69f2f89767p-40, -0x1.4ed5e276b35c3p-34, 0x1.519e17cc82153p-28, -0x1.61e0c93fbac5dp-22, // subdomain 3:6 (4.75..4.875)
                0x1.dcdb3ca17da7fp-63, -0x1.7bd4a2ad39c91p-58, 0x1.3f2984456d216p-52, -0x1.1f9f71e378d5fp-46, 0x1.132bda5ef6543p-40, -0x1.14bf39e3c08e1p-34, 0x1.2203d34693da7p-28, -0x1.3abed7bfaecc8p-22, // subdomain 3:7 (4.875..5.0)
                0x1.351712364b852p-63, -0x1.02fd48c91b50ap-58, 0x1.cd831bb97808fp-53, -0x1.b70e4542765fbp-47, 0x1.b83531bc78410p-41, -0x1.cd182aef322afp-35, 0x1.f5089e0082de2p-29, -0x1.190990a1bc197p-22, // subdomain 3:8 (5.0..5.125)
                0x1.966d1b8a14dcdp-64, -0x1.6862b8598316fp-59, 0x1.53cb45edff65ap-53, -0x1.53bbc8c619ffcp-47, 0x1.6394e053395cep-41, -0x1.82d3f9203923cp-35, 0x1.b2fefdac99aefp-29, -0x1.f7b1920eefcfdp-23, // subdomain 3:9 (5.125..5.25)
                0x1.107d541ef7742p-64, -0x1.ffb6d4b5dcf5bp-60, 0x1.fc8c953bca3acp-54, -0x1.0a0875b83c32ap-47, 0x1.21b70b3e10832p-41, -0x1.46953d36f265fp-35, 0x1.7b6be3a0ea4c7p-29, -0x1.c4e3f8c137de7p-23, // subdomain 3:10 (5.25..5.375)
                0x1.7548c6a941537p-65, -0x1.72262fda4797ep-60, 0x1.81ff8db0bac4bp-54, -0x1.a500ffe981901p-48, 0x1.dbbeec5dbbdfep-42, -0x1.154fa5ee79ef0p-35, 0x1.4c5e5fb3d1310p-29, -0x1.987bd82b3cdfep-23, // subdomain 3:11 (5.375..5.5)
                0x1.050ca2a94552ap-65, -0x1.10394eb2e8b95p-60, 0x1.289eb0812c74fp-54, -0x1.503336030b4ffp-48, 0x1.89555e8ea4202p-42, -0x1.d96f7fd431882p-36, 0x1.244f86399dd16p-29, -0x1.7181e0ab968dep-23, // subdomain 3:12 (5.5..5.625)
                0x1.74143110b0625p-66, -0x1.963c9b75150e1p-61, 0x1.ccd3b7e739ad0p-55, -0x1.0eb0bf70723eap-48, 0x1.4741991c3e8aap-42, -0x1.961a990f79fabp-36, 0x1.0208c29c72367p-29, -0x1.4f2a60c144be0p-23, // subdomain 3:13 (5.625..5.75)
                0x1.0d9a2d2074ca4p-66, -0x1.32f2589482c83p-61, 0x1.696886534a424p-55, -0x1.b721f3f18762fp-49, 0x1.11dc5986af828p-42, -0x1.5dec796368c1fp-36, 0x1.c9224fd266e7dp-30, -0x1.30cc748320582p-23, // subdomain 3:14 (5.75..5.875)
                0x1.8c602a94fd459p-67, -0x1.d4fefe1f385a5p-62, 0x1.1de03d383aa76p-55, -0x1.669c1e956ade9p-49, 0x1.cccff4baa94cap-43, -0x1.2ec942c985efdp-36, 0x1.963f405a37632p-30, -0x1.15db388b3aec8p-23, // subdomain 3:15 (5.875..6.0)
                0x1.2711ffdb4a8f1p-67, -0x1.69d4057a5a06dp-62, 0x1.c7c7bf24c0f5ep-56, -0x1.26ad0b160a680p-49, 0x1.859f8ce72f170p-43, -0x1.0707fa6f3703dp-36, 0x1.6a212279186f9p-30, -0x1.fbc0f91e2ee25p-24, // subdomain 3:16 (6.0..6.125)
                0x1.bc36a934ffabcp-68, -0x1.199ea9fa98b0fp-62, 0x1.6deac5e3cf1b1p-56, -0x1.e711ba2529641p-50, 0x1.4af60935aefb8p-43, -0x1.caaea429e6dc0p-37, 0x1.43bc02321d622p-30, -0x1.d0f12e9244534p-24, // subdomain 3:17 (6.125..6.25)
                0x1.51b2c97c50a84p-68, -0x1.b9e894dca62fcp-63, 0x1.27b24aafd852ep-56, -0x1.94b146db5da68p-50, 0x1.1a5a4e3a9a3b4p-43, -0x1.91531e611381ep-37, 0x1.22331eb4645e2p-30, -0x1.aa9e2e151ca59p-24, // subdomain 3:18 (6.25..6.375)
                0x1.0307f050492a3p-68, -0x1.5d450bd2fea2ep-63, 0x1.e0cdcd669e8aap-57, -0x1.51eda6c07a385p-50, 0x1.e3bbb27150459p-44, -0x1.604b0a52857f8p-37, 0x1.04cf95fe5a359p-30, -0x1.88383896e9af8p-24, // subdomain 3:19 (6.375..6.5)
                0x1.90a2e875905b2p-69, -0x1.15edef388b6c0p-63, 0x1.891abe578e018p-57, -0x1.1b800b5744608p-50, 0x1.9ff6d376de88bp-44, -0x1.3636776b82588p-37, 0x1.d5f2382dbf46ep-31, -0x1.694470c12305ap-24, // subdomain 3:20 (6.5..6.625)
                0x1.3829b2ef53b82p-69, -0x1.bd1dba503b581p-64, 0x1.4319b9a0b8fa6p-57, -0x1.ddc57dc9085dcp-51, 0x1.66fcd547d8100p-44, -0x1.11f68cc620309p-37, 0x1.a8607f6c7c714p-31, -0x1.4d59678537535p-24, // subdomain 3:21 (6.625..6.75)
                0x1.e9d9a7fe95cc8p-70, -0x1.66898efb1da74p-64, 0x1.0ae1545741b88p-57, -0x1.944180fed10b4p-51, 0x1.36e13a5b51781p-44, -0x1.e5406c44a6e89p-38, 0x1.8013bec6703ebp-31, -0x1.341c4c3a4e818p-24, // subdomain 3:22 (6.75..6.875)
                0x1.82d646d0cc2cfp-70, -0x1.2264ab3a86705p-64, 0x1.baf191782f12bp-58, -0x1.5764fe3889177p-51, 0x1.0e1898e4fdf6bp-44, -0x1.aee39c8fc52fap-38, 0x1.5c5713791bf46p-31, -0x1.1d3e9ec569b17p-24, // subdomain 3:23 (6.875..7.0)
                0x1.3357b95f54de3p-70, -0x1.d8d7999d10e87p-65, 0x1.713300bfcd1b8p-58, -0x1.24c7bdc7a21f8p-51, 0x1.d6c74c55ec499p-45, -0x1.7f9617fb3e50bp-38, 0x1.3c91299fa3d25p-31, -0x1.087c49cf1f1e4p-24, // subdomain 3:24 (7.0..7.125)
                0x1.eb270197418f3p-71, -0x1.82d7dc6366447p-65, 0x1.350412e2d36d1p-58, -0x1.f5039d11f5507p-52, 0x1.9b7e582d6e7a8p-45, -0x1.564d501ad5b27p-38, 0x1.203f89dc43fa1p-31, -0x1.eb342175b2f45p-25, // subdomain 3:25 (7.125..7.25)
                0x1.8a8d3dfd25a49p-71, -0x1.3df2cf951a2c8p-65, 0x1.03aa89a473e9dp-58, -0x1.ae1d11bbc725bp-52, 0x1.68af54e8019c2p-45, -0x1.322aac35ef66cp-38, 0x1.06f2c973796d2p-31, -0x1.c8c882d2e480dp-25, // subdomain 3:26 (7.25..7.375)
                0x1.3e8e6f46e8ed8p-71, -0x1.0677cc8c06b56p-65, 0x1.b60906e7561f0p-59, -0x1.726d445b114f0p-52, 0x1.3cff77cd80f5dp-45, -0x1.1273874a37620p-38, 0x1.e096dc5ae2461p-32, -0x1.a95b3a4ab309cp-25, // subdomain 3:27 (7.375..7.5)
                0x1.026e8e9d8571cp-71, -0x1.b3237f8fd7abap-66, 0x1.72c94efbb04d9p-59, -0x1.3ffef08a50177p-52, 0x1.1751281f678fep-45, -0x1.ed158cf88f2dcp-39, 0x1.b7eec386a30e9p-32, -0x1.8c9cdd53915e9p-25, // subdomain 3:28 (7.5..7.625)
                0x1.a537b68ead40bp-72, -0x1.6a20f7602c6a9p-66, 0x1.3aeeb1aa7ff63p-59, -0x1.153c7f912abb2p-52, 0x1.ed706199492a9p-46, -0x1.bbd75f5b82459p-39, 0x1.935fb11423ce6p-32, -0x1.72478a6fc7ff7p-25, // subdomain 3:29 (7.625..7.75)
                0x1.58c13c5bccc56p-72, -0x1.2e811e73a4a03p-66, 0x1.0c5d075ccf2d4p-59, -0x1.e1b876912a7b3p-53, 0x1.b4e0d4dfa4efap-46, -0x1.904b7901669b7p-39, 0x1.726f2a354c022p-32, -0x1.5a1d9c7dfec82p-25, // subdomain 3:30 (7.75..7.875)
                0x1.1b556a669182ap-72, -0x1.fb36248d1328ap-67, 0x1.cac923aeff3fep-60, -0x1.a3a1d58444652p-53, 0x1.83ac9e3b0e3b6p-46, -0x1.69b26afd3e05dp-39, 0x1.54b2eec01f2c7p-32, -0x1.43e890b4baf9dp-25, // subdomain 3:31 (7.875..8.0)
                0x1.a990e7813dc8ap-63, -0x1.8817855ffc640p-58, 0x1.6c64c32d5ae08p-52, -0x1.56ce5f34f4de7p-46, 0x1.455f2deee8080p-40, -0x1.37b55d6bb23dap-34, 0x1.2d53bce58eca8p-28, -0x1.25dbf63307042p-22, // subdomain 4:0 (8.0..8.25)
                0x1.25f2bb414192fp-63, -0x1.190b2e422a0f7p-58, 0x1.0edab9e938ce8p-52, -0x1.07fe5e5cf6e99p-46, 0x1.036dc6beec3bbp-40, -0x1.0127000d63e83p-34, 0x1.010baedb07df4p-28, -0x1.030cb50cd6989p-22, // subdomain 4:1 (8.25..8.5)
                0x1.9b92e62f3ded9p-64, -0x1.97b1f1031acbap-59, 0x1.96d5a00b1d2b0p-53, -0x1.9a3ea9d0507e1p-47, 0x1.a0da3928ed432p-41, -0x1.aafbd4401fa51p-35, 0x1.b8cd46512adbfp-29, -0x1.ca8f7985184d1p-23, // subdomain 4:2 (8.5..8.75)
                0x1.23bb04d62e661p-64, -0x1.2af7ae3811e9ep-59, 0x1.347b0aed180e5p-53, -0x1.4169704d31349p-47, 0x1.514674f547a8ap-41, -0x1.649662a3ec342p-35, 0x1.7bc7881017605p-29, -0x1.97637b1783a5ep-23, // subdomain 4:3 (8.75..9.0)
                0x1.a25810b89075ep-65, -0x1.baf3b5b1bd39fp-60, 0x1.d7fe2d21d4d94p-54, -0x1.fb88979bf3406p-48, 0x1.12b1eeb2056d5p-41, -0x1.2b7235efd0902p-35, 0x1.48adee6a84093p-29, -0x1.6b33edd62ce56p-23, // subdomain 4:4 (9.0..9.25)
                0x1.2f2c2868e17ecp-65, -0x1.4b42a75f56eebp-60, 0x1.6c19af1f2dcd1p-54, -0x1.939e9bfa9bfb7p-48, 0x1.c2395487a87a7p-42, -0x1.f9887f4b764a9p-36, 0x1.1da8a2aa7d162p-29, -0x1.44e1f4235a80fp-23, // subdomain 4:5 (9.25..9.5)
                0x1.bbcf5167c2466p-66, -0x1.f3dcf354c348fp-61, 0x1.1b1003a67b2a7p-54, -0x1.4328d484a0173p-48, 0x1.7319f988f8f26p-42, -0x1.accf022d96ee7p-36, 0x1.f283d7030cb3ap-30, -0x1.2383444cf89acp-23, // subdomain 4:6 (9.5..9.75)
                0x1.47e1317e0d906p-66, -0x1.7c466f7c8cea2p-61, 0x1.bb57e9f98df4ap-55, -0x1.0463472565df3p-48, 0x1.338e363b6726cp-42, -0x1.6d661b38d717cp-36, 0x1.b49f5f2016d80p-30, -0x1.06575e2db71eep-23, // subdomain 4:7 (9.75..10.0)
                0x1.e8b583f5220eap-67, -0x1.238e15a88e5a2p-61, 0x1.5d9439fae54a8p-55, -0x1.a622a5ab3c20cp-49, 0x1.00341c3e17b58p-42, -0x1.38b6ffd1fce2fp-36, 0x1.7fc52a7f4f4f9p-30, -0x1.d97e63caf4756p-24, // subdomain 4:8 (10.0..10.25)
                0x1.6f364e055e085p-67, -0x1.c2598f3224346p-62, 0x1.156e79584ba05p-55, -0x1.581d500d63023p-49, 0x1.aced5536a5f89p-43, -0x1.0cb916fefece4p-36, 0x1.527279dcc649fp-30, -0x1.ac6d3b2461509p-24, // subdomain 4:9 (10.25..10.5)
                0x1.160f463b461b4p-67, -0x1.5e394c5806f87p-62, 0x1.bb0acd11e988cp-56, -0x1.1a0378ffd2e70p-49, 0x1.68b2edbdff0b9p-43, -0x1.cfa1717a534dbp-37, 0x1.2b6d03cfdd70fp-30, -0x1.849ea1d8b333fp-24, // subdomain 4:10 (10.5..10.75)
                0x1.a82f641195639p-68, -0x1.1223dfa8bcab6p-62, 0x1.63cda14fe70fbp-56, -0x1.d0943ab4ba189p-50, 0x1.30a54c0ca5b86p-43, -0x1.916ca47bf857ep-37, 0x1.09b3f8095ea44p-30, -0x1.615985f6547c6p-24, // subdomain 4:11 (10.75..11.0)
                0x1.45ca3b02d0892p-68, -0x1.afd3bab52d4a2p-63, 0x1.1f4f0c35b828bp-56, -0x1.808081382be42p-50, 0x1.025e179dd0459p-43, -0x1.5cc7544bf9da9p-37, 0x1.d8e8c8268f36cp-31, -0x1.4201f68556961p-24, // subdomain 4:12 (11.0..11.25)
                0x1.f7b7d185da8f9p-69, -0x1.561ac24ee8960p-63, 0x1.d269e2b70c339p-57, -0x1.3fad8aeac48c1p-50, 0x1.b7f63a58fa79fp-44, -0x1.300b4de28f7d5p-37, 0x1.a6004d6679582p-31, -0x1.2613fe816fbe7p-24, // subdomain 4:13 (11.25..11.5)
                0x1.87d513fabc0d4p-69, -0x1.1088e4d2e6324p-63, 0x1.7c74f42e67ed9p-57, -0x1.0aef55e0f7b7ap-50, 0x1.77ff2e645c379p-44, -0x1.09e2dbb9ab9cbp-37, 0x1.798dac38716cap-31, -0x1.0d1f80d663ddap-24, // subdomain 4:14 (11.5..11.75)
                0x1.3299b48bfaf57p-69, -0x1.b486eea8a4bb1p-64, 0x1.37cd383c1109dp-57, -0x1.bfa0dc2e4f1dep-51, 0x1.427a3bb1f72e5p-44, -0x1.d270099438623p-38, 0x1.52a049b9be3b5p-31, -0x1.ed89bd3fd428ap-25, // subdomain 4:15 (11.75..12.0)
            };

            alignas(64) const double coeffs1[144 * 4] = {
                -0x1.9d5e2f887fe99p-15, 0x1.461380c17af85p-8, 0x1.ef90904c7eeeep-2, 0, // subdomain 0:0 (0.5..0.515625)
                -0x1.a0aac55ff7b75p-15, 0x1.391aec16a4582p-8, 0x1.f98d5b983b626p-2, 0, // subdomain 0:1 (0.515625..0.53125)
                -0x1.a34eb4d108597p-15, 0x1.2c0a96301c7a0p-8, 0x1.0190ffdb18ffap-1, 0, // subdomain 0:2 (0.53125..0.546875)
                -0x1.a54b75905c7b1p-15, 0x1.1ee7bea414b01p-8, 0x1.0626e9d159016p-1, 0, // subdomain 0:3 (0.546875..0.5625)
                -0x1.a6a2e09f9c0dap-15, 0x1.11b797c3f9012p-8, 0x1.0a882c11d80b0p-1, 0, // subdomain 0:4 (0.5625..0.578125)
                -0x1.a7572d47d55b1p-15, 0x1.047f439def882p-8, 0x1.0eb49ba94b015p-1, 0, // subdomain 0:5 (0.578125..0.59375)
                -0x1.a76aedcf732c1p-15, 0x1.ee87a22f32788p-9, 0x1.12ac220740fddp-1, 0, // subdomain 0:6 (0.59375..0.609375)
                -0x1.a6e10bf0336f4p-15, 0x1.d41472484c1d9p-9, 0x1.166ebcabe3f22p-1, 0, // subdomain 0:7 (0.609375..0.625)
                -0x1.a5bcc511c4ff5p-15, 0x1.b9aeb82f53c22p-9, 0x1.19fc7ccacb7eep-1, 0, // subdomain 0:8 (0.625..0.640625)
                -0x1.a401a64ddfab0p-15, 0x1.9f5ffe32f3cfdp-9, 0x1.1d5586e35a9eap-1, 0, // subdomain 0:9 (0.640625..0.65625)
                -0x1.a1b38842cf5acp-15, 0x1.8531936e65bbep-9, 0x1.207a124f27936p-1, 0, // subdomain 0:10 (0.65625..0.671875)
                -0x1.9ed68ab98a0d7p-15, 0x1.6b2c8733c03f0p-9, 0x1.236a68c6f3b80p-1, 0, // subdomain 0:11 (0.671875..0.6875)
                -0x1.9b6f1024805eap-15, 0x1.5159a4bd25507p-9, 0x1.2626e5debe65bp-1, 0, // subdomain 0:12 (0.6875..0.703125)
                -0x1.9781b8fc695b9p-15, 0x1.37c16f271a477p-9, 0x1.28aff67983042p-1, 0, // subdomain 0:13 (0.703125..0.71875)
                -0x1.93135f0054a9ap-15, 0x1.1e6c1db6001fcp-9, 0x1.2b061835369e2p-1, 0, // subdomain 0:14 (0.71875..0.734375)
                -0x1.8e29105e517edp-15, 0x1.056198684d07bp-9, 0x1.2d29d8cf9ce8fp-1, 0, // subdomain 0:15 (0.734375..0.75)
                -0x1.88c80ac9f5c35p-15, 0x1.d952e9ada710ep-10, 0x1.2f1bd5848f9cap-1, 0, // subdomain 0:16 (0.75..0.765625)
                -0x1.82f5b6860718dp-15, 0x1.a895e6c81e4a0p-10, 0x1.30dcba66553d8p-1, 0, // subdomain 0:17 (0.765625..0.78125)
                -0x1.7cb7a1667799fp-15, 0x1.7899f97635587p-10, 0x1.326d41b0a6140p-1, 0, // subdomain 0:18 (0.78125..0.796875)
                -0x1.761379cfd233fp-15, 0x1.496c3f3d43ff8p-10, 0x1.33ce3316ff061p-1, 0, // subdomain 0:19 (0.796875..0.8125)
                -0x1.6f0f09b916b9dp-15, 0x1.1b191e0fa05c0p-10, 0x1.3500630ee256fp-1, 0, // subdomain 0:20 (0.8125..0.828125)
                -0x1.67b031b4e48c3p-15, 0x1.db588253aef1cp-11, 0x1.3604b216a6047p-1, 0, // subdomain 0:21 (0.828125..0.84375)
                -0x1.5ffce406ac395p-15, 0x1.82612d1c18fa8p-11, 0x1.36dc0bf96e90ap-1, 0, // subdomain 0:22 (0.84375..0.859375)
                -0x1.57fb1fc8751b8p-15, 0x1.2b609a376898ep-11, 0x1.37876710f384fp-1, 0, // subdomain 0:23 (0.859375..0.875)
                -0x1.4fb0ec25941e9p-15, 0x1.acd34d0552fc5p-12, 0x1.3807c385a8ea6p-1, 0, // subdomain 0:24 (0.875..0.890625)
                -0x1.472453ae7cc93p-15, 0x1.071b591598bb4p-12, 0x1.385e2a8de655bp-1, 0, // subdomain 0:25 (0.890625..0.90625)
                -0x1.3e5b5fc99dd8ap-15, 0x1.96e4222f45e07p-14, 0x1.388badaca1013p-1, 0, // subdomain 0:26 (0.90625..0.921875)
                -0x1.355c1444ff8c7p-15, -0x1.b9b7c13bd75dap-15, 0x1.389165f04ad0ep-1, 0, // subdomain 0:27 (0.921875..0.9375)
                -0x1.2c2c6b0c1b92ep-15, -0x1.9f35f85e34ca3p-13, 0x1.3870733264021p-1, 0, // subdomain 0:28 (0.9375..0.953125)
                -0x1.22d2500525cd3p-15, -0x1.635c505586adbp-12, 0x1.3829fb5848c3ep-1, 0, // subdomain 0:29 (0.953125..0.96875)
                -0x1.19539d18ba3cdp-15, -0x1.f26733025ac65p-12, 0x1.37bf2995bf073p-1, 0, // subdomain 0:30 (0.96875..0.984375)
                -0x1.0fb616669fcb8p-15, -0x1.3e55652364b90p-11, 0x1.37312db1c3897p-1, 0, // subdomain 0:31 (0.984375..1.0)
                -0x1.011c5cf8a63b1p-13, -0x1.a17e5468f1efbp-10, 0x1.361ce6aff8f20p-1, 0, // subdomain 1:0 (1.0..1.03125)
                -0x1.dad97f327156bp-14, -0x1.0e91436da8346p-9, 0x1.343cc43eebb11p-1, 0, // subdomain 1:1 (1.03125..1.0625)
                -0x1.b3479ed963454p-14, -0x1.47736f828d6f2p-9, 0x1.31e5ec7ae1c19p-1, 0, // subdomain 1:2 (1.0625..1.09375)
                -0x1.8bd2d546498e1p-14, -0x1.7b64635f21d42p-9, 0x1.2f2242328728ep-1, 0, // subdomain 1:3 (1.09375..1.125)
                -0x1.64c56720344d7p-14, -0x1.aa6c756c73d75p-9, 0x1.2bfba10c072f2p-1, 0, // subdomain 1:4 (1.125..1.15625)
                -0x1.3e63a6c82a173p-14, -0x1.d49ce6a50e315p-9, 0x1.287bcaf670341p-1, 0, // subdomain 1:5 (1.15625..1.1875)
                -0x1.18eb9610d846dp-14, -0x1.fa0f1e1cc0cfdp-9, 0x1.24ac57171221ep-1, 0, // subdomain 1:6 (1.1875..1.21875)
                -0x1.e9294776a07b7p-15, -0x1.0d71ed3ecf96fp-8, 0x1.2096a24a6c3e0p-1, 0, // subdomain 1:7 (1.21875..1.25)
                -0x1.a31f07936f638p-15, -0x1.1ba12e727ee57p-8, 0x1.1c43c14950261p-1, 0, // subdomain 1:8 (1.25..1.28125)
                -0x1.600c42ab1f475p-15, -0x1.27abc8337f4ddp-8, 0x1.17bc747c2c5e4p-1, 0, // subdomain 1:9 (1.28125..1.3125)
                -0x1.20375082c804dp-15, -0x1.31aa98a7b74acp-8, 0x1.13091d8011b8ep-1, 0, // subdomain 1:10 (1.3125..1.34375)
                -0x1.c7b2935886a5fp-16, -0x1.39b87b7e885b1p-8, 0x1.0e31b65b1182ap-1, 0, // subdomain 1:11 (1.34375..1.375)
                -0x1.563cdc8d9290ep-16, -0x1.3ff1e17962469p-8, 0x1.093dca5818a73p-1, 0, // subdomain 1:12 (1.375..1.40625)
                -0x1.d89abe7ee37c8p-17, -0x1.44746bb90b20cp-8, 0x1.0434707d7c211p-1, 0, // subdomain 1:13 (1.40625..1.4375)
                -0x1.14195f590513fp-17, -0x1.475e8bd343d56p-8, 0x1.fe388f1a1f348p-2, 0, // subdomain 1:14 (1.4375..1.46875)
                -0x1.7c6ca88d3d7e0p-19, -0x1.48cf298d05830p-8, 0x1.f3f6e6f1be1a6p-2, 0, // subdomain 1:15 (1.46875..1.5)
                0x1.198a66e7bcb7fp-19, -0x1.48e54eedb36a1p-8, 0x1.e9af386906071p-2, 0, // subdomain 1:16 (1.5..1.53125)
                0x1.b90eb89ae9483p-18, -0x1.47bfdb380fd70p-8, 0x1.df6bdb9109ec1p-2, 0, // subdomain 1:17 (1.53125..1.5625)
                0x1.639b8bbe32409p-17, -0x1.457d3d2f3dbefp-8, 0x1.d5363316c7294p-2, 0, // subdomain 1:18 (1.5625..1.59375)
                0x1.dc04c31376bb1p-17, -0x1.423b34e93d5c9p-8, 0x1.cb16b0c4efe99p-2, 0, // subdomain 1:19 (1.59375..1.625)
                0x1.2320f8ebc1916p-16, -0x1.3e169d4da9bd4p-8, 0x1.c114dbd87a161p-2, 0, // subdomain 1:20 (1.625..1.65625)
                0x1.5173c5ee79890p-16, -0x1.392b3d3f79e75p-8, 0x1.b73758ea816c7p-2, 0, // subdomain 1:21 (1.65625..1.6875)
                0x1.794de76ec47eap-16, -0x1.3393a05383ccep-8, 0x1.ad83f3328b595p-2, 0, // subdomain 1:22 (1.6875..1.71875)
                0x1.9b0994fe0a0b4p-16, -0x1.2d68f6dbbb6f3p-8, 0x1.a3ffa6e596920p-2, 0, // subdomain 1:23 (1.71875..1.75)
                0x1.b70668bb4bf10p-16, -0x1.26c2fcf8d0025p-8, 0x1.9aaeac787c47dp-2, 0, // subdomain 1:24 (1.75..1.78125)
                0x1.cda7c38a18b26p-16, -0x1.1fb7e850f03a3p-8, 0x1.9194848ee97bfp-2, 0, // subdomain 1:25 (1.78125..1.8125)
                0x1.df534f0e88935p-16, -0x1.185c5bfb2695cp-8, 0x1.88b404638064fp-2, 0, // subdomain 1:26 (1.8125..1.84375)
                0x1.ec6f9ffd6cfa1p-16, -0x1.10c36222d4aebp-8, 0x1.800f627968365p-2, 0, // subdomain 1:27 (1.84375..1.875)
                0x1.f562fa821b29bp-16, -0x1.08fe6ade376bfp-8, 0x1.77a8436b99fdbp-2, 0, // subdomain 1:28 (1.875..1.90625)
                0x1.fa9239bada958p-16, -0x1.011d4fad5dfb3p-8, 0x1.6f7fc6b279f96p-2, 0, // subdomain 1:29 (1.90625..1.9375)
                0x1.fc5fda9c76fdap-16, -0x1.f25cb628c00a9p-9, 0x1.6796933cadcaep-2, 0, // subdomain 1:30 (1.9375..1.96875)
                0x1.fb2b29eea441bp-16, -0x1.e27ca78742ba4p-9, 0x1.5fece3bc8b246p-2, 0, // subdomain 1:31 (1.96875..2.0)
                0x1.f47e9306cfe2ap-14, -0x1.cad949ac74a7ep-8, 0x1.54e50a117281cp-2, 0, // subdomain 2:0 (2.0..2.0625)
                0x1.e43f494374dbdp-14, -0x1.ac098c3390d17p-8, 0x1.470a2c751e645p-2, 0, // subdomain 2:1 (2.0625..2.125)
                0x1.cdc8ee1532684p-14, -0x1.8e7276a1783f1p-8, 0x1.3a212c60c909fp-2, 0, // subdomain 2:2 (2.125..2.1875)
                0x1.b32b5c23639a0p-14, -0x1.7266851f59e5ap-8, 0x1.2e1ee4a5eed45p-2, 0, // subdomain 2:3 (2.1875..2.25)
                0x1.961e21d988c32p-14, -0x1.5819fa15bdf64p-8, 0x1.22f618ca6a7c2p-2, 0, // subdomain 2:4 (2.25..2.3125)
                0x1.7803ad4d576e4p-14, -0x1.3fa84ebd0404bp-8, 0x1.189850f22ad9fp-2, 0, // subdomain 2:5 (2.3125..2.375)
                0x1.59efbc877bdd5p-14, -0x1.2919552fb6a1cp-8, 0x1.0ef68b5b88c36p-2, 0, // subdomain 2:6 (2.375..2.4375)
                0x1.3cafecda1c526p-14, -0x1.1465e110bc5bdp-8, 0x1.0601c69a216cdp-2, 0, // subdomain 2:7 (2.4375..2.5)
                0x1.20d56b8d013dcp-14, -0x1.017bdcfd44971p-8, 0x1.fb56d19a94ff8p-3, 0, // subdomain 2:8 (2.5..2.5625)
                0x1.06bef8e1ca385p-14, -0x1.e08382db6b806p-9, 0x1.ebcb1143aa4e6p-3, 0, // subdomain 2:9 (2.5625..2.625)
                0x1.dd453e1d68e0fp-15, -0x1.c132da26b969bp-9, 0x1.dd463a3414083p-3, 0, // subdomain 2:10 (2.625..2.6875)
                0x1.b12d5a7fb419ep-15, -0x1.a4c4ce5c33027p-9, 0x1.cfb031e56850ep-3, 0, // subdomain 2:11 (2.6875..2.75)
                0x1.893360481682ep-15, -0x1.8af734012e8b7p-9, 0x1.c2f2ec3f9f193p-3, 0, // subdomain 2:12 (2.75..2.8125)
                0x1.653288afc6418p-15, -0x1.73892da2a56e4p-9, 0x1.b6fa6ab751a86p-3, 0, // subdomain 2:13 (2.8125..2.875)
                0x1.44ef74ee8e31ap-15, -0x1.5e3cec77eca99p-9, 0x1.abb4aa62a05a1p-3, 0, // subdomain 2:14 (2.875..2.9375)
                0x1.282131a02139fp-15, -0x1.4ad8cede3aa4dp-9, 0x1.a11186a9d6be2p-3, 0, // subdomain 2:15 (2.9375..3.0)
                0x1.0e7828826d6eap-15, -0x1.3927ff7c36d9ep-9, 0x1.97029517a7c68p-3, 0, // subdomain 2:16 (3.0..3.0625)
                0x1.ef4686dec632ep-16, -0x1.28fab41f5fd34p-9, 0x1.8d7afdc4af746p-3, 0, // subdomain 2:17 (3.0625..3.125)
                0x1.c6a71453f9529p-16, -0x1.1a2626f6ae0f7p-9, 0x1.846f52f2cb861p-3, 0, // subdomain 2:18 (3.125..3.1875)
                0x1.a27d08f887226p-16, -0x1.0c84602d90bdep-9, 0x1.7bd569aa440d4p-3, 0, // subdomain 2:19 (3.1875..3.25)
                0x1.823f4aac9f477p-16, -0x1.ffe7c2e7f2a18p-10, 0x1.73a4348cb11f8p-3, 0, // subdomain 2:20 (3.25..3.3125)
                0x1.6571f7dafda2fp-16, -0x1.e8ae81dfd0636p-10, 0x1.6bd3a195b87ebp-3, 0, // subdomain 2:21 (3.3125..3.375)
                0x1.4ba6bd7ba850bp-16, -0x1.d3297b199f0c5p-10, 0x1.645c7b245c666p-3, 0, // subdomain 2:22 (3.375..3.4375)
                0x1.347c75beefa95p-16, -0x1.bf2ba75553c05p-10, 0x1.5d384c6420abap-3, 0, // subdomain 2:23 (3.4375..3.5)
                0x1.1f9e52478bdecp-16, -0x1.ac8dada0e1199p-10, 0x1.566148ff153ccp-3, 0, // subdomain 2:24 (3.5..3.5625)
                0x1.0cc2c4b5165c9p-16, -0x1.9b2d2660dbb3bp-10, 0x1.4fd237e46297ep-3, 0, // subdomain 2:25 (3.5625..3.625)
                0x1.f75494e0d5ac3p-17, -0x1.8aebf0df81b2fp-10, 0x1.498660de3f32cp-3, 0, // subdomain 2:26 (3.625..3.6875)
                0x1.d83c6aaf7b17ap-17, -0x1.7baf9c53a1fc2p-10, 0x1.43797caac8e6fp-3, 0, // subdomain 2:27 (3.6875..3.75)
                0x1.bbdf03f2d74f8p-17, -0x1.6d60e40188bd5p-10, 0x1.3da7a74a45ce0p-3, 0, // subdomain 2:28 (3.75..3.8125)
                0x1.a1eb8a8bd3738p-17, -0x1.5feb3d46092f8p-10, 0x1.380d54388f9c7p-3, 0, // subdomain 2:29 (3.8125..3.875)
                0x1.8a1c6d64e826fp-17, -0x1.533c75e6e80a6p-10, 0x1.32a7444d3017bp-3, 0, // subdomain 2:30 (3.875..3.9375)
                0x1.7435b0b1fe346p-17, -0x1.474460d1822e7p-10, 0x1.2d727d04bf46dp-3, 0, // subdomain 2:31 (3.9375..4.0)
                0x1.567fc34962e0bp-15, -0x1.3687a8af074d2p-9, 0x1.25f9cb1d2d3c5p-3, 0, // subdomain 3:0 (4.0..4.125)
                0x1.33e2a62f841a5p-15, -0x1.223b8ef41aa26p-9, 0x1.1c982f4493a06p-3, 0, // subdomain 3:1 (4.125..4.25)
                0x1.16042e9c508cep-15, -0x1.0ff210d378a36p-9, 0x1.13d0b73a39076p-3, 0, // subdomain 3:2 (4.25..4.375)
                0x1.f81eb863de7a3p-16, -0x1.fecc662f6ada6p-10, 0x1.0b946aeb2e515p-3, 0, // subdomain 3:3 (4.375..4.5)
                0x1.caba5ed77f9d6p-16, -0x1.e0bd8fb61ea2cp-10, 0x1.03d648ff217bfp-3, 0, // subdomain 3:4 (4.5..4.625)
                0x1.a2d1ae896575ap-16, -0x1.c557e5396eaffp-10, 0x1.f915e5b6b75edp-4, 0, // subdomain 3:5 (4.625..4.75)
                0x1.7f8eaf2717f23p-16, -0x1.ac4a97dd0f4c6p-10, 0x1.eb50d3c384e8cp-4, 0, // subdomain 3:6 (4.75..4.875)
                0x1.604317e7c6611p-16, -0x1.9550ed48207c6p-10, 0x1.de4bb3668bb37p-4, 0, // subdomain 3:7 (4.875..5.0)
                0x1.445f912ca8000p-16, -0x1.80300e2d5afd7p-10, 0x1.d1f6d8de9ca8cp-4, 0, // subdomain 3:8 (5.0..5.125)
                0x1.2b6d34dae4a8bp-16, -0x1.6cb54db32d4cfp-10, 0x1.c6444d76b0754p-4, 0, // subdomain 3:9 (5.125..5.25)
                0x1.1508a4f63fa06p-16, -0x1.5ab4cb6dfde27p-10, 0x1.bb2793deb136cp-4, 0, // subdomain 3:10 (5.25..5.375)
                0x1.00de46f225bb8p-16, -0x1.4a085a4f24705p-10, 0x1.b0957654732d0p-4, 0, // subdomain 3:11 (5.375..5.5)
                0x1.dd4eaba768d78p-17, -0x1.3a8e9cd25bd72p-10, 0x1.a683dcb74c707p-4, 0, // subdomain 3:12 (5.5..5.625)
                0x1.bc4f293f3ee62p-17, -0x1.2c2a4ac9e77ebp-10, 0x1.9ce9a90d1b30fp-4, 0, // subdomain 3:13 (5.625..5.75)
                0x1.9e56f45164287p-17, -0x1.1ec1980cfd413p-10, 0x1.93be995054e18p-4, 0, // subdomain 3:14 (5.75..5.875)
                0x1.830d73e907b39p-17, -0x1.123db5606cf42p-10, 0x1.8afb2d9bd97cdp-4, 0, // subdomain 3:15 (5.875..6.0)
                0x1.6a262aa14dc70p-17, -0x1.068a667787d78p-10, 0x1.829891f81f746p-4, 0, // subdomain 3:16 (6.0..6.125)
                0x1.535ec9e1fed6fp-17, -0x1.f72b521abb009p-11, 0x1.7a908b319cabcp-4, 0, // subdomain 3:17 (6.125..6.25)
                0x1.3e7d9e716233dp-17, -0x1.e29ed3e1f1a1fp-11, 0x1.72dd663cbd983p-4, 0, // subdomain 3:18 (6.25..6.375)
                0x1.2b5044298a4eap-17, -0x1.cf528b0bc7a24p-11, 0x1.6b79e9c22e522p-4, 0, // subdomain 3:19 (6.375..6.5)
                0x1.19aa92bd00acep-17, -0x1.bd2ca37507776p-11, 0x1.64614980181a8p-4, 0, // subdomain 3:20 (6.5..6.625)
                0x1.0965b88df84cep-17, -0x1.ac15dfb31c88ap-11, 0x1.5d8f1b3b4a17dp-4, 0, // subdomain 3:21 (6.625..6.75)
                0x1.f4bef5e6e698bp-18, -0x1.9bf949d3e4affp-11, 0x1.56ff4d06d5237p-4, 0, // subdomain 3:22 (6.75..6.875)
                0x1.d8f32a26d17d9p-18, -0x1.8cc3ef1c68b69p-11, 0x1.50ae1cb207886p-4, 0, // subdomain 3:23 (6.875..7.0)
                0x1.bf325fb7b8a55p-18, -0x1.7e64a50a7b8c8p-11, 0x1.4a981034622f8p-4, 0, // subdomain 3:24 (7.0..7.125)
                0x1.a74cdd71fa7bep-18, -0x1.70cbd62ba4d99p-11, 0x1.44b9eef56e319p-4, 0, // subdomain 3:25 (7.125..7.25)
                0x1.91181afd8d3d4p-18, -0x1.63eb559ad1bf0p-11, 0x1.3f10bbd38d83ap-4, 0, // subdomain 3:26 (7.25..7.375)
                0x1.7c6e194e9b592p-18, -0x1.57b63829543aap-11, 0x1.3999afd131e12p-4, 0, // subdomain 3:27 (7.375..7.5)
                0x1.692cd31cca4cap-18, -0x1.4c20b26141d44p-11, 0x1.345235537d117p-4, 0, // subdomain 3:28 (7.5..7.625)
                0x1.5735c17ec90bep-18, -0x1.411ffab136817p-11, 0x1.2f37e3e04f3c3p-4, 0, // subdomain 3:29 (7.625..7.75)
                0x1.466d7181606b3p-18, -0x1.36aa2f2d64300p-11, 0x1.2a487c4c4e553p-4, 0, // subdomain 3:30 (7.75..7.875)
                0x1.36bb281d50870p-18, -0x1.2cb63e68d9842p-11, 0x1.2581e54b928efp-4, 0, // subdomain 3:31 (7.875..8.0)
                0x1.2108c824a6046p-16, -0x1.1ea9c38006723p-10, 0x1.1ea0476d73f3ep-4, 0, // subdomain 4:0 (8.0..8.25)
                0x1.0728f01297f67p-16, -0x1.0d6c8f152a0e4p-10, 0x1.15f10209dda72p-4, 0, // subdomain 4:1 (8.25..8.5)
                0x1.e09c44d0f43c3p-17, -0x1.fb691d5828beep-11, 0x1.0dc57147e40e7p-4, 0, // subdomain 4:2 (8.5..8.75)
                0x1.b814204d94df9p-17, -0x1.dea9fef384120p-11, 0x1.0612232a78184p-4, 0, // subdomain 4:3 (8.75..9.0)
                0x1.9401cbd4f6c41p-17, -0x1.c44ed4ebe46a2p-11, 0x1.fd99e3adbcae0p-5, 0, // subdomain 4:4 (9.0..9.25)
                0x1.73cb850419eb7p-17, -0x1.ac153427fa40dp-11, 0x1.efd9ab0ac8dc4p-5, 0, // subdomain 4:5 (9.25..9.5)
                0x1.56efded08ebc3p-17, -0x1.95c384791de6cp-11, 0x1.e2d37be675bebp-5, 0, // subdomain 4:6 (9.5..9.75)
                0x1.3d015ac3ed8eep-17, -0x1.81279fd4c6f19p-11, 0x1.d678e3df95dbfp-5, 0, // subdomain 4:7 (9.75..10.0)
                0x1.25a2e67c715d4p-17, -0x1.6e15b08a1e9cap-11, 0x1.cabce7d4dcf8cp-5, 0, // subdomain 4:8 (10.0..10.25)
                0x1.10850a7223ef3p-17, -0x1.5c6741da99764p-11, 0x1.bf93d53db675ap-5, 0, // subdomain 4:9 (10.25..10.5)
                0x1.fac746d54efa6p-18, -0x1.4bfa79241c58ap-11, 0x1.b4f31a564b07cp-5, 0, // subdomain 4:10 (10.5..10.75)
                0x1.d808125a8fcf6p-18, -0x1.3cb16fe7aac17p-11, 0x1.aad123fa36ff6p-5, 0, // subdomain 4:11 (10.75..11.0)
                0x1.b86717b0d21bfp-18, -0x1.2e71a8956213bp-11, 0x1.a12540405abf2p-5, 0, // subdomain 4:12 (11.0..11.25)
                0x1.9b8c6caa70f12p-18, -0x1.21239945a7fc3p-11, 0x1.97e78515b5f8bp-5, 0, // subdomain 4:13 (11.25..11.5)
                0x1.812b9086a74eep-18, -0x1.14b2487f5263bp-11, 0x1.8f10ba37f66a4p-5, 0, // subdomain 4:14 (11.5..11.75)
                0x1.6901b94e7875dp-18, -0x1.090af8ec5da93p-11, 0x1.869a460c3e163p-5, 0, // subdomain 4:15 (11.75..12.0)
            };
            // clang-format on
        } // namespace imwofx

        namespace dimwofx
        {
            // clang-format off
            alignas(64) const double coeffs0[72 * 8] = {
                0x1.d6462e2ba17ffp-60, -0x1.e551ecc8e41d5p-54, 0x1.200524e8cc7dfp-48, 0x1.42860b02ca18ep-43, -0x1.aeba5b90b3e38p-36, 0x1.5a510b24c63a4p-30, -0x1.8033fd13af310p-26, -0x1.2fa0a60c1045fp-20, // subdomain 0:8 (2.5..2.5625)
                0x1.13114e4b83dc3p-59, -0x1.94a679c730767p-54, 0x1.450e696adbdd9p-49, 0x1.b24a3d7f053adp-43, -0x1.84c48493bf721p-36, 0x1.0d4c5528a4128p-30, -0x1.8077fda21f85ap-27, -0x1.532738588c138p-20, // subdomain 0:9 (2.5625..2.625)
                0x1.1a1f5587536bep-59, -0x1.3cda8b1b976eep-54, 0x1.e80f3acee564bp-51, 0x1.e93a43e6152d0p-43, -0x1.51ccd2d6fbfe8p-36, 0x1.922d8541b1d89p-31, -0x1.6de56ec19e70dp-29, -0x1.61950c8c3936fp-20, // subdomain 0:10 (2.625..2.6875)
                0x1.08e1dd619abbdp-59, -0x1.cdbdf414a8b88p-55, -0x1.fbc8c6c068490p-53, 0x1.f396276fc53bbp-43, -0x1.1b5fc18f0e8a7p-36, 0x1.1dac45662b654p-31, 0x1.e9a8e29236925p-29, -0x1.603c6dccb9f41p-20, // subdomain 0:11 (2.6875..2.75)
                0x1.cf7c8b678f395p-60, -0x1.31eebcf13ed05p-55, -0x1.163f3444ee019p-50, 0x1.dd640f04e4ae3p-43, -0x1.cb66ee7b785dcp-37, 0x1.7b187d7ff1ae6p-32, 0x1.0ddf45e5efb0ep-27, -0x1.53aa0a6e8f466p-20, // subdomain 0:12 (2.75..2.8125)
                0x1.7c49ac62bb322p-60, -0x1.5ba7c6705ab6ep-56, -0x1.9bc1b7579b41fp-50, 0x1.b18e35f6d16a3p-43, -0x1.6776a44117416p-37, 0x1.c3b7d88e58106p-33, 0x1.6b6b76d3a0580p-27, -0x1.3f9fe1b903ac3p-20, // subdomain 0:13 (2.8125..2.875)
                0x1.243acfae71890p-60, -0x1.1301684ff3935p-57, -0x1.de7a230b10edfp-50, 0x1.7959fe2236e8ap-43, -0x1.0e9e6b0a6a8d6p-37, 0x1.b06a2b837a7d0p-34, 0x1.9ebbc97ee35d3p-27, -0x1.271d9f5c04d29p-20, // subdomain 0:14 (2.875..2.9375)
                0x1.a06d7646aa17ep-61, 0x1.28038500c1355p-60, -0x1.edc6c684d0f76p-50, 0x1.3c2aca1eed75dp-43, -0x1.85761eb698252p-38, 0x1.544fea4da4cd0p-36, 0x1.b234fa3133f8ep-27, -0x1.0c71f3dfd3bc3p-20, // subdomain 0:15 (2.9375..3.0)
                0x1.0c22d9bf9351bp-61, 0x1.f369d0e90e5ebp-58, -0x1.d859e10e122d3p-50, 0x1.feff9983718c4p-44, -0x1.0881f52a7466fp-38, -0x1.3cffe57386c24p-35, 0x1.aeba979a802e1p-27, -0x1.e2a24189f8e29p-21, // subdomain 0:16 (3.0..3.0625)
                0x1.23b36de8e9431p-62, 0x1.79a7fdf33fad1p-57, -0x1.ab322436eb133p-50, 0x1.8e36862326043p-44, -0x1.4a89b2b21fe29p-39, -0x1.3e1f7ed9440a6p-34, 0x1.9ba05dffe720dp-27, -0x1.addb6169c1d80p-21, // subdomain 0:17 (3.0625..3.125)
                0x1.9c2df56d2ed27p-64, 0x1.b5e7c10f42ad5p-57, -0x1.710bff6b36e74p-50, 0x1.2a85bad981531p-44, -0x1.6550cc30494aap-40, -0x1.9c3bda035fd31p-34, 0x1.7eb9f0b6b2351p-27, -0x1.7c222fcc53fa6p-21, // subdomain 0:18 (3.125..3.1875)
                -0x1.00736717780c5p-65, 0x1.bfb083476e4bep-57, -0x1.3233cff52d26ap-50, 0x1.ac2ec68898469p-45, -0x1.0c633d125dea5p-41, -0x1.c9288d1126c95p-34, 0x1.5c7e5d4ccdd68p-27, -0x1.4e655b41da62dp-21, // subdomain 0:19 (3.1875..3.25)
                -0x1.dfebef67406a2p-64, 0x1.a6d327d687c6ap-57, -0x1.e93d46b1233f5p-51, 0x1.229cf6b2afe5fp-45, 0x1.5eb606f85a679p-44, -0x1.d2be70c1a347dp-34, 0x1.3836a04d79904p-27, -0x1.2518163548f32p-21, // subdomain 0:20 (3.25..3.3125)
                -0x1.526454b34adb9p-63, 0x1.78c3614cad1b8p-57, -0x1.7860dae307f7dp-51, 0x1.6d520852bf7acp-46, 0x1.f1acf47d536b0p-42, -0x1.c4497efa0d8a7p-34, 0x1.143059731cf7cp-27, -0x1.0054b17ea05e3p-21, // subdomain 0:21 (3.3125..3.375)
                -0x1.798580f92c1cbp-63, 0x1.402f6d8672986p-57, -0x1.16343851ce7d3p-51, 0x1.949222587b0bap-47, 0x1.73362d586b6a2p-41, -0x1.a6aa5c8a01616p-34, 0x1.e3df847c3882ep-28, -0x1.bff1c30a8f0f1p-22, // subdomain 0:22 (3.375..3.4375)
                -0x1.76872b33935a7p-63, 0x1.04fe1e24c075ap-57, -0x1.88fa9b73c0822p-52, 0x1.50ffcca1f78a6p-48, 0x1.b066d2b11b4e2p-41, -0x1.8097d71f9ef79p-34, 0x1.a4bb51147606ep-28, -0x1.8777fb5611b84p-22, // subdomain 0:23 (3.4375..3.5)
                -0x1.585a4c84ca972p-63, 0x1.99219c36b0428p-58, -0x1.064caf5eab90dp-52, 0x1.7940894f86cb9p-53, 0x1.c24901afc434cp-41, -0x1.56f09281839cap-34, 0x1.6bdf25325581bp-28, -0x1.567fb4e306304p-22, // subdomain 0:24 (3.5..3.5625)
                -0x1.2b12eeb1839a2p-63, 0x1.3456a9cb014b9p-58, -0x1.4395e8b2a6698p-53, -0x1.8c2dc03f91e36p-49, 0x1.b72e627ca8072p-41, -0x1.2d112a958ac45p-34, 0x1.399289643b6fap-28, -0x1.2c3a103533be1p-22, // subdomain 0:25 (3.5625..3.625)
                -0x1.ef92fe8568521p-64, 0x1.bda94d87a06abp-59, -0x1.5e245cf9586a5p-54, -0x1.40ea789e83893p-48, 0x1.9a1f2757e1ab6p-41, -0x1.0526bfd9c9cacp-34, 0x1.0daf4aa7b0a54p-28, -0x1.07d696db74138p-22, // subdomain 0:26 (3.625..3.6875)
                -0x1.89a65a85b0f26p-64, 0x1.32dbe72eee9fap-59, -0x1.187964e8dcdc5p-55, -0x1.7cb90d082a35cp-48, 0x1.7345e00d01f43p-41, -0x1.c0f20056ac56ep-35, 0x1.cf95f66435ad1p-29, -0x1.d11c77eff4ea7p-23, // subdomain 0:27 (3.6875..3.75)
                -0x1.2c585e810ec10p-64, 0x1.8d6d56aa26d3ap-60, 0x1.726d08bea9f2ap-64, -0x1.8d0bfbe25e9afp-48, 0x1.486b2ee65cacbp-41, -0x1.7f525e74bc644p-35, 0x1.8eab71e6d6edep-29, -0x1.9b53b9d90a346p-23, // subdomain 0:28 (3.75..3.8125)
                -0x1.b7da8923bfa4ep-65, 0x1.d72157158d6a9p-61, 0x1.6003790b8bc9cp-56, -0x1.8122236fb2fa8p-48, 0x1.1d719f4246e74p-41, -0x1.45caf639f89d6p-35, 0x1.576009a2e3b25p-29, -0x1.6d0afa1a2670cp-23, // subdomain 0:29 (3.8125..3.875)
                -0x1.33e47cdca4c1bp-65, 0x1.dde949acaa38dp-62, 0x1.119c4a1a4ed6cp-55, -0x1.64690df842d87p-48, 0x1.e98f3981e5abdp-42, -0x1.1422387913e35p-35, 0x1.2872efca1cd19p-29, -0x1.45227a6f2ea99p-23, // subdomain 0:30 (3.875..3.9375)
                -0x1.987f429072a52p-66, 0x1.413d9a8f48f20p-63, 0x1.3cbd2c4056064p-55, -0x1.3f16c7beff50ep-48, 0x1.9f9431c52c89cp-42, -0x1.d389db37ffd66p-36, 0x1.00b400d2cf2dcp-29, -0x1.22a1b0e700c09p-23, // subdomain 0:31 (3.9375..4.0)
                -0x1.7f590b1bad16ap-57, -0x1.bc133ae52cf55p-55, 0x1.3e54bb7068f67p-47, -0x1.02b9a4f8b190fp-41, 0x1.40c71d91320fbp-36, -0x1.6c74336b9bc30p-31, 0x1.a06d8c64c388cp-26, -0x1.ee77ea69f4f7ap-21, // subdomain 1:0 (4.0..4.125)
                -0x1.0e6f931bcad34p-59, -0x1.de708ecc4c53dp-54, 0x1.07b809a9d434ap-47, -0x1.729d939e6ffe7p-42, 0x1.c04c07f0b073fp-37, -0x1.0749319fe221dp-31, 0x1.3f4946aab560cp-26, -0x1.9329572650c7dp-21, // subdomain 1:1 (4.125..4.25)
                0x1.d5fd56cda5022p-60, -0x1.d60747971546ap-54, 0x1.85c9e97b5f932p-48, -0x1.000994db0884bp-42, 0x1.3879aef090ef5p-37, -0x1.817e48468e851p-32, 0x1.f162df0b630aap-27, -0x1.4c9f4e4e964b3p-21, // subdomain 1:2 (4.25..4.375)
                0x1.636df67dbc6bdp-59, -0x1.736d2c43760e1p-54, 0x1.0e74e632a7781p-48, -0x1.5c199b7c18bd8p-43, 0x1.b6f6ec912e2abp-38, -0x1.1ededdc760768p-32, 0x1.893c4142274ddp-27, -0x1.15475c12ee96cp-21, // subdomain 1:3 (4.375..4.5)
                0x1.3d80a4df10bf6p-59, -0x1.081ccc5a16e90p-54, 0x1.6ad726f1de85bp-49, -0x1.d895f901a2982p-44, 0x1.38c03d4b1dfbfp-38, -0x1.b251d810be32cp-33, 0x1.3b18664605611p-27, -0x1.d278116fda26dp-22, // subdomain 1:4 (4.5..4.625)
                0x1.d97b60a59d854p-60, -0x1.628235f99d32fp-55, 0x1.e030669365353p-50, -0x1.43b86a1b055d6p-44, 0x1.c58b711ee992bp-39, -0x1.4e2c8aa39aa6ep-33, 0x1.fef3707a68404p-28, -0x1.8b7916cab9312p-22, // subdomain 1:5 (4.625..4.75)
                0x1.4240fb66f80b7p-60, -0x1.cdfe4fdbc67a4p-56, 0x1.3e1a22da9f7dfp-50, -0x1.c283b07ccaf6ap-45, 0x1.4eee40d996584p-39, -0x1.04ecf2a77a431p-33, 0x1.a28b5b136b4f8p-28, -0x1.519e17cca6287p-22, // subdomain 1:6 (4.75..4.875)
                0x1.a193f8543aa49p-61, -0x1.2a1a12571d5f3p-56, 0x1.aa4c9a6e93c12p-51, -0x1.3f7401d222859p-45, 0x1.f7572354b2efbp-40, -0x1.9cc1c17230f00p-34, 0x1.59ef085c1207ep-28, -0x1.2203d346ab214p-22, // subdomain 1:7 (4.875..5.0)
                0x1.0867d4873887cp-61, -0x1.8273094e1f8e4p-57, 0x1.22b931f608c86p-51, -0x1.cde3acce8d61bp-46, 0x1.802c8e5907a3bp-40, -0x1.4a27e1576ea10p-34, 0x1.202f1ad51ae45p-28, -0x1.f5089e00a10afp-23, // subdomain 1:8 (5.0..5.125)
                0x1.4d658ebb1feefp-62, -0x1.fc242e927c93dp-58, 0x1.94a0955d6cd20p-52, -0x1.540ac32421c94p-46, 0x1.29445add7535ap-40, -0x1.0aafa5a3cd300p-34, 0x1.e388f767c8bf8p-29, -0x1.b2fefdacad85cp-23, // subdomain 1:9 (5.125..5.25)
                0x1.a82c2578c8ad7p-63, -0x1.54adce37d3e2ep-58, 0x1.1f537ea859368p-52, -0x1.fce1b7aa85f67p-47, 0x1.d18edc3e5c7cap-41, -0x1.b2928d5f3540dp-35, 0x1.983a8c845e69fp-29, -0x1.7b6be3a0f7999p-23, // subdomain 1:10 (5.25..5.375)
                0x1.12917c679a6b3p-63, -0x1.d2b0295ed17d5p-59, 0x1.9fc0eaa8d9b9ep-53, -0x1.8239de2302124p-47, 0x1.7060e9230ff69p-41, -0x1.64cf2ee203f18p-35, 0x1.5aa38f69e4414p-29, -0x1.4c5e5fb3da4d7p-23, // subdomain 1:11 (5.375..5.5)
                0x1.6b310e01eeb67p-64, -0x1.465d125ed3f3ap-59, 0x1.31d0019389dc5p-53, -0x1.28c778a57466ap-47, 0x1.262cd55ab561ap-41, -0x1.2700053ec7f09p-35, 0x1.27e5afe47c753p-29, -0x1.244f8639a4309p-23, // subdomain 1:12 (5.5..5.625)
                0x1.eb40693a96f6ep-65, -0x1.d12a38c98a197p-60, 0x1.c86c1050af402p-54, -0x1.cd0dd8b6db6aap-48, 0x1.d9b55742e5436p-42, -0x1.eae263480847fp-36, 0x1.fba13f5329cf6p-30, -0x1.0208c29c76c10p-23, // subdomain 1:13 (5.625..5.75)
                0x1.5346c3fa9bb6ep-65, -0x1.510bd5940281fp-60, 0x1.58e794755ee70p-54, -0x1.6992a4e1a6bbep-48, 0x1.803dbb24a9b68p-42, -0x1.9aca848fc8610p-36, 0x1.b56797bc22b80p-30, -0x1.c9224fd26d7c9p-24, // subdomain 1:14 (5.75..5.875)
                0x1.dd7f0457be63dp-66, -0x1.ef8716dcfc812p-61, 0x1.078580d79a88fp-54, -0x1.1dff33b4a72b6p-48, 0x1.39c89ec44220cp-42, -0x1.599bf646e42ddp-36, 0x1.7a7b937bd0be0p-30, -0x1.963f405a3c39ap-24, // subdomain 1:15 (5.875..6.0)
                0x1.558d09380be82p-66, -0x1.70e0abe2f06b0p-61, 0x1.96a4c2a1eb0b0p-55, -0x1.c7f5d88d5171dp-49, 0x1.01d76c90d5c52p-42, -0x1.2437a8bb5ed13p-36, 0x1.48c9f90af48bep-30, -0x1.6a2122791c099p-24, // subdomain 1:16 (6.0..6.125)
                0x1.ef9ce723286f9p-67, -0x1.15a93e9affea9p-61, 0x1.3c85c3072886cp-55, -0x1.6e0d792ed1843p-49, 0x1.aa2f8708f2391p-43, -0x1.f0710c642d245p-37, 0x1.1ead269a24847p-30, -0x1.43bc023220182p-24, // subdomain 1:17 (6.125..6.25)
                0x1.6c16446e24671p-67, -0x1.a629805d0b877p-62, 0x1.f0b4e96e27ca4p-56, -0x1.27ccabf69bf9dp-49, 0x1.621b210de5915p-43, -0x1.a7877442ead33p-37, 0x1.f5a7e5f947177p-31, -0x1.22331eb4666dbp-24, // subdomain 1:18 (6.25..6.375)
                0x1.0e6eebb8dc744p-67, -0x1.43d11b995ea96p-62, 0x1.8899ef682eb8ep-56, -0x1.e0f645a6520f7p-50, 0x1.27aff42d3aafdp-43, -0x1.6accc50085876p-37, 0x1.b85dcce71a074p-31, -0x1.04cf95fe5bca4p-24, // subdomain 1:19 (6.375..6.5)
                0x1.95bd9bbf27f5dp-68, -0x1.f4d6110a14374p-63, 0x1.386cdb19f43ccp-56, -0x1.893a0a596d2a6p-50, 0x1.f0201740214cdp-44, -0x1.37f91df4d8520p-37, 0x1.83c41546594bfp-31, -0x1.d5f2382dc1b8dp-25, // subdomain 1:20 (6.5..6.625)
                0x1.331ff35f927f2p-68, -0x1.863bc79cbbf72p-63, 0x1.f4625711d94bap-57, -0x1.43321c5b72236p-50, 0x1.a20cd0a384bd3p-44, -0x1.0d3d9f75dc0ddp-37, 0x1.56742ff7a0f19p-31, -0x1.a8607f6c7e58fp-25, // subdomain 1:21 (6.625..6.75)
                0x1.d4ce509d6f476p-69, -0x1.322db5b12d383p-63, 0x1.93122b27b40a7p-57, -0x1.0af4767391787p-50, 0x1.61b952d664ed7p-44, -0x1.d251d6c0148a7p-38, 0x1.2f4843aae2c05p-31, -0x1.8013bec671bd6p-25, // subdomain 1:22 (6.75..6.875)
                0x1.6889317454f57p-69, -0x1.e394556931a79p-64, 0x1.46796e58826a5p-57, -0x1.bb0fc99fc3c4bp-51, 0x1.2c785ff4a1d52p-44, -0x1.9524e4b8d681cp-38, 0x1.0d4e41d9d6f5dp-31, -0x1.5c5713791d229p-25, // subdomain 1:23 (6.875..7.0)
                0x1.173d26bd03bfap-69, -0x1.80340e47fb78cp-64, 0x1.09ce0a6d88ab9p-57, -0x1.714b032934a78p-51, 0x1.002ec73a8aa20p-44, -0x1.611578c264d5dp-38, 0x1.df7b9dfa07434p-32, -0x1.3c91299fa4c27p-25, // subdomain 1:24 (7.0..7.125)
                0x1.b3697853f5117p-70, -0x1.32fd3e22bba64p-64, 0x1.b2ef6dd6b2dc4p-58, -0x1.35174212f5e74p-51, 0x1.b6632b434894fp-45, -0x1.349ec1bd5b55ep-38, 0x1.abe0a42185f39p-32, -0x1.203f89dc44b9ep-25, // subdomain 1:25 (7.125..7.25)
                0x1.559258a251df4p-70, -0x1.ed37ff68241a2p-65, 0x1.657c46e8a33e8p-58, -0x1.03b9f2e71c37ap-51, 0x1.785970f31ad35p-45, -0x1.0e837f5d18c10p-38, 0x1.7eb5574367322p-32, -0x1.06f2c9737a074p-25, // subdomain 1:26 (7.25..7.375)
                0x1.0d84f73f2ebb7p-70, -0x1.8e37c78f31521p-65, 0x1.271d0b5fdd926p-58, -0x1.b621e9a2b0e24p-52, 0x1.441f9cf11c09ep-45, -0x1.db7f33319b59bp-39, 0x1.5710691cc2074p-32, -0x1.e096dc5ae33efp-26, // subdomain 1:27 (7.375..7.5)
                0x1.abac308fc4fb9p-71, -0x1.430ea4ce33fe4p-65, 0x1.e945b8298304cp-59, -0x1.72dd7f5042852p-52, 0x1.17ff135ea74afp-45, -0x1.a2f9bbc51dcdep-39, 0x1.342d781b56f2ap-32, -0x1.b7eec386a3d87p-26, // subdomain 1:28 (7.5..7.625)
                0x1.55149067708cdp-71, -0x1.07464a11d3aabp-65, 0x1.973047255fd24p-59, -0x1.3aff259a5a14cp-52, 0x1.e529e0ac505b1p-46, -0x1.721448dc965bcp-39, 0x1.15669b992f651p-32, -0x1.935fb1142472ep-26, // subdomain 1:29 (7.625..7.75)
                0x1.115cb3f00ac4cp-71, -0x1.aef6fc4efb393p-66, 0x1.5426ef33d990bp-59, -0x1.0c6a7eb91b316p-52, 0x1.a58168e49283bp-46, -0x1.47a89f61092dcp-39, 0x1.f45e5741bd035p-33, -0x1.726f2a354c88dp-26, // subdomain 1:30 (7.75..7.875)
                0x1.b83ad4b3acb46p-72, -0x1.622f0f35810b8p-66, 0x1.1d2c600652f81p-59, -0x1.cadf460f85a11p-53, 0x1.6f2d9bc01b225p-46, -0x1.22c1767230875p-39, 0x1.c41f05bc8aea3p-33, -0x1.54b2eec01f9b1p-26, // subdomain 1:31 (7.875..8.0)
                0x1.40fd2dfa32de0p-62, -0x1.0a06b3c4683ebp-57, 0x1.b853ab1f26d09p-52, -0x1.6ca73e77da328p-46, 0x1.2bf49e13df6c9p-40, -0x1.e80ebf7247780p-35, 0x1.85a2b4c624d69p-29, -0x1.2d53bce5a390cp-23, // subdomain 2:0 (8.0..8.25)
                0x1.aaeedc6269661p-63, -0x1.6f7ef13accb23p-58, 0x1.3ba86026ae52cp-52, -0x1.0f08a5b1d46bcp-46, 0x1.cdfd33765e244p-41, -0x1.8524a65a105a3p-35, 0x1.4170c0106bc66p-29, -0x1.010baedb1638cp-23, // subdomain 2:1 (8.25..8.5)
                0x1.2043793e9fc41p-63, -0x1.0145ea68dbd9ep-58, 0x1.c9f5a8953707dp-53, -0x1.9715ec35176c8p-47, 0x1.66f6de42ae818p-41, -0x1.38a3a83b95c6dp-35, 0x1.0add64a7dd037p-29, -0x1.b8cd46513ef38p-24, // subdomain 2:2 (8.5..8.75)
                0x1.8abcb25898f50p-64, -0x1.6cb720f723131p-59, 0x1.4fdc6850e3e2bp-53, -0x1.34a89e4df8849p-47, 0x1.193c48e317dccp-41, -0x1.f9e9abb2d93c6p-36, 0x1.bdbbfb4c9c432p-30, -0x1.7bc78810259e4p-24, // subdomain 2:3 (8.75..9.0)
                0x1.11c5eb063c647p-64, -0x1.057ffeb93c590p-59, 0x1.f1a89f3074d03p-54, -0x1.d83f886baf6c4p-48, 0x1.bc178dd8623afp-42, -0x1.9c0ae35cccc25p-36, 0x1.764ec36b90b12p-30, -0x1.48adee6a8e3f6p-24, // subdomain 2:4 (9.0..9.25)
                0x1.80517fb9151e4p-65, -0x1.7b035d6a80520p-60, 0x1.7433fa8149676p-54, -0x1.6c490c5792e73p-48, 0x1.612aceee156c9p-42, -0x1.51aafd746daf3p-36, 0x1.3bf54f8f056cep-30, -0x1.1da8a2aa847cap-24, // subdomain 2:5 (9.25..9.5)
                0x1.10c3d4e0939d7p-65, -0x1.1569f18847a0bp-60, 0x1.18d7d22aeb72ap-54, -0x1.1b32aeb11c958p-48, 0x1.1ac3be87aeb11p-42, -0x1.165379bab83edp-36, 0x1.0c01615c646c0p-30, -0x1.f283d7031788fp-25, // subdomain 2:6 (9.5..9.75)
                0x1.873899b0bd40bp-66, -0x1.99e523b2ca7b9p-61, 0x1.ab56186aaf19ep-55, -0x1.bb8b238dbc1b2p-49, 0x1.c7adc311e645ep-43, -0x1.cd554f3f40326p-37, 0x1.c8bfa206e7b41p-31, -0x1.b49f5f201ed8fp-25, // subdomain 2:7 (9.75..10.0)
                0x1.1b4c7d891331fp-66, -0x1.3179a3bc0bdb1p-61, 0x1.47a81e9b067ccp-55, -0x1.5dba6707bd950p-49, 0x1.715e55b6b6e30p-43, -0x1.804e28cc4bfa0p-37, 0x1.86e4bfc661328p-31, -0x1.7fc52a7f55469p-25, // subdomain 2:8 (10.0..10.25)
                0x1.9e0ddfcd1be6ap-67, -0x1.cb0f8599ddaf0p-62, 0x1.fa248a09df928p-56, -0x1.158b28c8bc4d4p-49, 0x1.2d19a984fad38p-43, -0x1.41b1febbcb63bp-37, 0x1.4fe75cbeaad76p-31, -0x1.527279dccac55p-25, // subdomain 2:9 (10.25..10.5)
                0x1.312fceeb223c1p-67, -0x1.5b9b712db591bp-62, 0x1.89a1f50121cb1p-56, -0x1.bb363e4de20bap-50, 0x1.ed8618de98465p-44, -0x1.0e86316a6d89bp-37, 0x1.21c4e6ec6591dp-31, -0x1.2b6d03cfe0d5cp-25, // subdomain 2:10 (10.5..10.75)
                0x1.c58acb66dbc90p-68, -0x1.0923a9446639fp-62, 0x1.342222ed65e2dp-56, -0x1.63eec42e3dbd3p-50, 0x1.9681b72c36456p-44, -0x1.c8f7f0b70b85fp-38, 0x1.f5c7cd9ae0e34p-32, -0x1.09b3f809613afp-25, // subdomain 2:11 (10.75..11.0)
                0x1.5397f5abcb48dp-68, -0x1.97459b439b20ap-63, 0x1.e56509a3cb7f2p-57, -0x1.1f687f668a07dp-50, 0x1.507073ea7f81ap-44, -0x1.838d22617fae1p-38, 0x1.b3f9295ee8301p-32, -0x1.d8e8c8269330bp-26, // subdomain 2:12 (11.0..11.25)
                0x1.001eca640cc12p-68, -0x1.3ad95f9238022p-63, 0x1.808ecb5189590p-57, -0x1.d2913c2e37a96p-51, 0x1.17b7dbb37e64ep-44, -0x1.49f8aaf426bbcp-38, 0x1.7c0e215b27328p-32, -0x1.a6004d667c6b2p-26, // subdomain 2:13 (11.25..11.5)
                0x1.84fd7628cdeabp-69, -0x1.e9d3f693ef662p-64, 0x1.325dc754089dcp-57, -0x1.7c93902a1757ap-51, 0x1.d322d98d21a2ap-45, -0x1.19ff622a92ac8p-38, 0x1.4c5b92a80d471p-32, -0x1.798dac3873d0dp-26, // subdomain 2:14 (11.5..11.75)
                0x1.295395723c946p-69, -0x1.7f474f42253c4p-64, 0x1.eabbbaafeb32fp-58, -0x1.37e52bbf5c432p-51, 0x1.87acc32717828p-45, -0x1.e3b7588f76847p-39, 0x1.238605fc9c2dcp-32, -0x1.52a049b9c01a5p-26, // subdomain 2:15 (11.75..12.0)
            };

            alignas(64) const double coeffs1[72 * 4] = {
                0x1.b048fc84c3ec1p-14, -0x1.20d56b8d013dap-8, 0x1.017bdcfd44971p-3, 0, // subdomain 0:8 (2.5..2.5625)
                0x1.91f7a55c72031p-14, -0x1.06bef8e1ca383p-8, 0x1.e08382db6b806p-4, 0, // subdomain 0:9 (2.5625..2.625)
                0x1.715a38b0abc21p-14, -0x1.dd453e1d68e0cp-9, 0x1.c132da26b969bp-4, 0, // subdomain 0:10 (2.625..2.6875)
                0x1.5029b152cba18p-14, -0x1.b12d5a7fb419cp-9, 0x1.a4c4ce5c33027p-4, 0, // subdomain 0:11 (2.6875..2.75)
                0x1.2fa87bc621853p-14, -0x1.893360481682cp-9, 0x1.8af734012e8b7p-4, 0, // subdomain 0:12 (2.75..2.8125)
                0x1.10b5596e90edfp-14, -0x1.653288afc6417p-9, 0x1.73892da2a56e4p-4, 0, // subdomain 0:13 (2.8125..2.875)
                0x1.e7bc25cae7da8p-15, -0x1.44ef74ee8e31ap-9, 0x1.5e3cec77eca99p-4, 0, // subdomain 0:14 (2.875..2.9375)
                0x1.b2e1db2ddab54p-15, -0x1.282131a02139fp-9, 0x1.4ad8cede3aa4dp-4, 0, // subdomain 0:15 (2.9375..3.0)
                0x1.831877a5dc07bp-15, -0x1.0e7828826d6ebp-9, 0x1.3927ff7c36d9ep-4, 0, // subdomain 0:16 (3.0..3.0625)
                0x1.585760c3b5e0cp-15, -0x1.ef4686dec632fp-10, 0x1.28fab41f5fd34p-4, 0, // subdomain 0:17 (3.0625..3.125)
                0x1.3266bc4e52af0p-15, -0x1.c6a71453f952ap-10, 0x1.1a2626f6ae0f7p-4, 0, // subdomain 0:18 (3.125..3.1875)
                0x1.10f0f429afff6p-15, -0x1.a27d08f887227p-10, 0x1.0c84602d90bdep-4, 0, // subdomain 0:19 (3.1875..3.25)
                0x1.e7204e6d0abf2p-16, -0x1.823f4aac9f478p-10, 0x1.ffe7c2e7f2a18p-5, 0, // subdomain 0:20 (3.25..3.3125)
                0x1.b3b022093ea48p-16, -0x1.6571f7dafda30p-10, 0x1.e8ae81dfd0636p-5, 0, // subdomain 0:21 (3.3125..3.375)
                0x1.86b9fe73ee749p-16, -0x1.4ba6bd7ba850cp-10, 0x1.d3297b199f0c5p-5, 0, // subdomain 0:22 (3.375..3.4375)
                0x1.5f70d3874fa56p-16, -0x1.347c75beefa96p-10, 0x1.bf2ba75553c05p-5, 0, // subdomain 0:23 (3.4375..3.5)
                0x1.3d176e95c76f0p-16, -0x1.1f9e52478bdedp-10, 0x1.ac8dada0e1199p-5, 0, // subdomain 0:24 (3.5..3.5625)
                0x1.1f034c5221eadp-16, -0x1.0cc2c4b5165cap-10, 0x1.9b2d2660dbb3bp-5, 0, // subdomain 0:25 (3.5625..3.625)
                0x1.049d7cd80671ep-16, -0x1.f75494e0d5ac3p-11, 0x1.8aebf0df81b2fp-5, 0, // subdomain 0:26 (3.625..3.6875)
                0x1.dac476dc139eap-17, -0x1.d83c6aaf7b17ap-11, 0x1.7baf9c53a1fc2p-5, 0, // subdomain 0:27 (3.6875..3.75)
                0x1.b1bf6dd0365e7p-17, -0x1.bbdf03f2d74f8p-11, 0x1.6d60e40188bd5p-5, 0, // subdomain 0:28 (3.75..3.8125)
                0x1.8d68cee2a00c2p-17, -0x1.a1eb8a8bd3738p-11, 0x1.5feb3d46092f8p-5, 0, // subdomain 0:29 (3.8125..3.875)
                0x1.6d1a676d155d1p-17, -0x1.8a1c6d64e826fp-11, 0x1.533c75e6e80a6p-5, 0, // subdomain 0:30 (3.875..3.9375)
                0x1.504723f8c04bep-17, -0x1.7435b0b1fe346p-11, 0x1.474460d1822e7p-5, 0, // subdomain 0:31 (3.9375..4.0)
                0x1.2a90cd3c08f1ep-15, -0x1.567fc34962e07p-10, 0x1.3687a8af074d2p-5, 0, // subdomain 1:0 (4.0..4.125)
                0x1.00adb82d1d43bp-15, -0x1.33e2a62f8419cp-10, 0x1.223b8ef41aa26p-5, 0, // subdomain 1:1 (4.125..4.25)
                0x1.bca716d0e33fap-16, -0x1.16042e9c508c4p-10, 0x1.0ff210d378a36p-5, 0, // subdomain 1:2 (4.25..4.375)
                0x1.83ad72faffa3ap-16, -0x1.f81eb863de795p-11, 0x1.fecc662f6ada6p-6, 0, // subdomain 1:3 (4.375..4.5)
                0x1.53f823bd77620p-16, -0x1.caba5ed77f9ccp-11, 0x1.e0bd8fb61ea2cp-6, 0, // subdomain 1:4 (4.5..4.625)
                0x1.2bae9a7ca1acfp-16, -0x1.a2d1ae8965753p-11, 0x1.c557e5396eaffp-6, 0, // subdomain 1:5 (4.625..4.75)
                0x1.096896efcc164p-16, -0x1.7f8eaf2717f1ep-11, 0x1.ac4a97dd0f4c6p-6, 0, // subdomain 1:6 (4.75..4.875)
                0x1.d81e439f86353p-17, -0x1.604317e7c660ep-11, 0x1.9550ed48207c6p-6, 0, // subdomain 1:7 (4.875..5.0)
                0x1.a58e58f29a27bp-17, -0x1.445f912ca7ffep-11, 0x1.80300e2d5afd7p-6, 0, // subdomain 1:8 (5.0..5.125)
                0x1.79c52d8b33dcdp-17, -0x1.2b6d34dae4a8ap-11, 0x1.6cb54db32d4cfp-6, 0, // subdomain 1:9 (5.125..5.25)
                0x1.53aafa90e9e78p-17, -0x1.1508a4f63fa05p-11, 0x1.5ab4cb6dfde27p-6, 0, // subdomain 1:10 (5.25..5.375)
                0x1.325ce2206da85p-17, -0x1.00de46f225bb7p-11, 0x1.4a085a4f24705p-6, 0, // subdomain 1:11 (5.375..5.5)
                0x1.15216880b0eabp-17, -0x1.dd4eaba768d77p-12, 0x1.3a8e9cd25bd72p-6, 0, // subdomain 1:12 (5.5..5.625)
                0x1.f6bf9121e71d5p-18, -0x1.bc4f293f3ee62p-12, 0x1.2c2a4ac9e77ebp-6, 0, // subdomain 1:13 (5.625..5.75)
                0x1.c932aec4b0847p-18, -0x1.9e56f45164286p-12, 0x1.1ec1980cfd413p-6, 0, // subdomain 1:14 (5.75..5.875)
                0x1.a0c8d4d0d862fp-18, -0x1.830d73e907b38p-12, 0x1.123db5606cf42p-6, 0, // subdomain 1:15 (5.875..6.0)
                0x1.7cd0bad6a329ep-18, -0x1.6a262aa14dc70p-12, 0x1.068a667787d78p-6, 0, // subdomain 1:16 (6.0..6.125)
                0x1.5cb4e2edb33e9p-18, -0x1.535ec9e1fed6fp-12, 0x1.f72b521abb009p-7, 0, // subdomain 1:17 (6.125..6.25)
                0x1.3ff6a28fd57c3p-18, -0x1.3e7d9e716233dp-12, 0x1.e29ed3e1f1a1fp-7, 0, // subdomain 1:18 (6.25..6.375)
                0x1.262a2a712f43bp-18, -0x1.2b5044298a4eap-12, 0x1.cf528b0bc7a24p-7, 0, // subdomain 1:19 (6.375..6.5)
                0x1.0ef35490da444p-18, -0x1.19aa92bd00acep-12, 0x1.bd2ca37507776p-7, 0, // subdomain 1:20 (6.5..6.625)
                0x1.f4061b47d2fd0p-19, -0x1.0965b88df84cep-12, 0x1.ac15dfb31c88ap-7, 0, // subdomain 1:21 (6.625..6.75)
                0x1.ce2a725775c24p-19, -0x1.f4bef5e6e698bp-13, 0x1.9bf949d3e4affp-7, 0, // subdomain 1:22 (6.75..6.875)
                0x1.abddee281e8a2p-19, -0x1.d8f32a26d17d9p-13, 0x1.8cc3ef1c68b69p-7, 0, // subdomain 1:23 (6.875..7.0)
                0x1.8cba6eb6aead7p-19, -0x1.bf325fb7b8a55p-13, 0x1.7e64a50a7b8c8p-7, 0, // subdomain 1:24 (7.0..7.125)
                0x1.7067191846374p-19, -0x1.a74cdd71fa7bep-13, 0x1.70cbd62ba4d99p-7, 0, // subdomain 1:25 (7.125..7.25)
                0x1.5696621e2b60ap-19, -0x1.91181afd8d3d4p-13, 0x1.63eb559ad1bf0p-7, 0, // subdomain 1:26 (7.25..7.375)
                0x1.3f046bb806475p-19, -0x1.7c6e194e9b592p-13, 0x1.57b63829543aap-7, 0, // subdomain 1:27 (7.375..7.5)
                0x1.2975a5fead06fp-19, -0x1.692cd31cca4cap-13, 0x1.4c20b26141d44p-7, 0, // subdomain 1:28 (7.5..7.625)
                0x1.15b5a7d3d5ff9p-19, -0x1.5735c17ec90bep-13, 0x1.411ffab136817p-7, 0, // subdomain 1:29 (7.625..7.75)
                0x1.0396355e7f162p-19, -0x1.466d7181606b3p-13, 0x1.36aa2f2d64300p-7, 0, // subdomain 1:30 (7.75..7.875)
                0x1.e5dcd90f1876cp-20, -0x1.36bb281d50870p-13, 0x1.2cb63e68d9842p-7, 0, // subdomain 1:31 (7.875..8.0)
                0x1.b8c9f14c8a881p-18, -0x1.2108c824a6044p-12, 0x1.1ea9c38006723p-7, 0, // subdomain 2:0 (8.0..8.25)
                0x1.84930f9341e61p-18, -0x1.0728f01297f65p-12, 0x1.0d6c8f152a0e4p-7, 0, // subdomain 2:1 (8.25..8.5)
                0x1.57eb9b23d23aap-18, -0x1.e09c44d0f43c0p-13, 0x1.fb691d5828beep-8, 0, // subdomain 2:2 (8.5..8.75)
                0x1.318a9c51a2bd0p-18, -0x1.b814204d94df7p-13, 0x1.dea9fef384120p-8, 0, // subdomain 2:3 (8.75..9.0)
                0x1.1066f260a1ac7p-18, -0x1.9401cbd4f6c3fp-13, 0x1.c44ed4ebe46a2p-8, 0, // subdomain 2:4 (9.0..9.25)
                0x1.e752ee3507c1fp-19, -0x1.73cb850419eb6p-13, 0x1.ac153427fa40dp-8, 0, // subdomain 2:5 (9.25..9.5)
                0x1.b544e67374e88p-19, -0x1.56efded08ebc3p-13, 0x1.95c384791de6cp-8, 0, // subdomain 2:6 (9.5..9.75)
                0x1.89830d4492aeap-19, -0x1.3d015ac3ed8edp-13, 0x1.81279fd4c6f19p-8, 0, // subdomain 2:7 (9.75..10.0)
                0x1.631ecad837584p-19, -0x1.25a2e67c715d4p-13, 0x1.6e15b08a1e9cap-8, 0, // subdomain 2:8 (10.0..10.25)
                0x1.4151ec5b48fc9p-19, -0x1.10850a7223ef3p-13, 0x1.5c6741da99764p-8, 0, // subdomain 2:9 (10.25..10.5)
                0x1.2376f96286671p-19, -0x1.fac746d54efa6p-14, 0x1.4bfa79241c58ap-8, 0, // subdomain 2:10 (10.5..10.75)
                0x1.09032478bf5d6p-19, -0x1.d808125a8fcf5p-14, 0x1.3cb16fe7aac17p-8, 0, // subdomain 2:11 (10.75..11.0)
                0x1.e302f1c801e13p-20, -0x1.b86717b0d21bfp-14, 0x1.2e71a8956213bp-8, 0, // subdomain 2:12 (11.0..11.25)
                0x1.b91dfdc2279dbp-20, -0x1.9b8c6caa70f12p-14, 0x1.21239945a7fc3p-8, 0, // subdomain 2:13 (11.25..11.5)
                0x1.93af414195cc9p-20, -0x1.812b9086a74eep-14, 0x1.14b2487f5263bp-8, 0, // subdomain 2:14 (11.5..11.75)
                0x1.72274defdf1e9p-20, -0x1.6901b94e7875cp-14, 0x1.090af8ec5da93p-8, 0, // subdomain 2:15 (11.75..12.0)
            };
            // clang-format on
        } // namespace dimwofx

        namespace erfcxtab
        {
            // clang-format off
            alignas(64) const double coeffs0[208 * 8] = {
                -0x1.6aa3b0e9e58ddp-87, 0x1.976d16e1f2025p-78, -0x1.b1d46d20c5e6dp-68, 0x1.b2f6a130d628fp-58, -0x1.973633d72fea1p-48, 0x1.601229af0bae5p-38, -0x1.14db078bd8f1ap-28, 0x1.8337f8244cb90p-19, // subdomain 0:0 (0.125..0.12890625)
                -0x1.64af62846d0a8p-87, 0x1.911aa5c1c017cp-78, -0x1.ab83683de3dcap-68, 0x1.ad134e2a96552p-58, -0x1.92262a770a25ap-48, 0x1.5c1e79e73ce39p-38, -0x1.121ed904461b3p-28, 0x1.7ffd834d4c102p-19, // subdomain 0:1 (0.12890625..0.1328125)
                -0x1.5ed5668be744dp-87, 0x1.8ae2c49b8925fp-78, -0x1.a54b77d0c0ee9p-68, 0x1.a745eab58d0f9p-58, -0x1.8d27aa12ab6bdp-48, 0x1.58375c3dee077p-38, -0x1.0f6a8542d5fb6p-28, 0x1.7ccb37336fdddp-19, // subdomain 0:2 (0.1328125..0.13671875)
                -0x1.5915433137e8fp-87, 0x1.84c4fe0f4fd06p-78, -0x1.9f2c32850dc2bp-68, 0x1.a18e1fc2eddf6p-58, -0x1.883a715e52756p-48, 0x1.545ca52e7ec32p-38, -0x1.0cbdf34eeb677p-28, 0x1.79a0fc6bf067fp-19, // subdomain 0:3 (0.13671875..0.140625)
                -0x1.536e80f112452p-87, 0x1.7ec0dedc0bdc7p-78, -0x1.992530d7b8f50p-68, 0x1.9beb97b162f6ep-58, -0x1.835e401144976p-48, 0x1.508e29d647765p-38, -0x1.0a190a865010dp-28, 0x1.767ebbd66de0bp-19, // subdomain 0:4 (0.140625..0.14453125)
                -0x1.4de0aa8874f76p-87, 0x1.78d5f5d56f0b1p-78, -0x1.93360d0e87cdcp-68, 0x1.965dfe46bf9a6p-58, -0x1.7e92d6e18edb0p-48, 0x1.4ccbbff216f94p-38, -0x1.077bb29bf3178p-28, 0x1.73645e9bef153p-19, // subdomain 0:5 (0.14453125..0.1484375)
                -0x1.486b4ce9617aap-87, 0x1.7303d3d9dd3eep-78, -0x1.8d5e632fd77d5p-68, 0x1.90e500a9cfa92p-58, -0x1.79d7f77fd9d21p-48, 0x1.49153ddbbae4bp-38, -0x1.04e5d396ac96ap-28, 0x1.7051ce2de3d95p-19, // subdomain 0:6 (0.1484375..0.15234375)
                -0x1.430df72fce521p-87, 0x1.6d4a0bc896e38p-78, -0x1.879dd0fa86126p-68, 0x1.8b804d5c43bc9p-58, -0x1.752d64934fd27p-48, 0x1.456a7a879230bp-38, -0x1.025755d006161p-28, 0x1.6d46f4452b211p-19, // subdomain 0:7 (0.15234375..0.15625)
                -0x1.3dc83a96d2972p-87, 0x1.67a8327814923p-78, -0x1.81f3f5de0263dp-68, 0x1.862f9434b9537p-58, -0x1.7092e1b595422p-48, 0x1.41cb4d8229fe1p-38, -0x1.ffa043e60f991p-29, 0x1.6a43bae11cb9fp-19, // subdomain 0:8 (0.15625..0.16015625)
                -0x1.3899aa6e09b65p-87, 0x1.621ddeac92f31p-78, -0x1.7c6072f2822b3p-68, 0x1.80f28658de829p-58, -0x1.6c08336ed2a26p-48, 0x1.3e378eede45bfp-38, -0x1.faa041f6153c8p-29, 0x1.67480c469694ap-19, // subdomain 0:9 (0.16015625..0.1640625)
                -0x1.3381dc0f2e2e0p-87, 0x1.5caaa90ecde18p-78, -0x1.76e2eaf15d90ap-68, 0x1.7bc8d637b0a61p-58, -0x1.678d1f31d0060p-48, 0x1.3aaf1780a8df0p-38, -0x1.f5ae78651d87ep-29, 0x1.6453d2ff0d929p-19, // subdomain 0:10 (0.1640625..0.16796875)
                -0x1.2e8066d3ea273p-87, 0x1.574e2c22e9deep-78, -0x1.717b022d8f6d5p-68, 0x1.76b23783d5924p-58, -0x1.63216b5821a4dp-48, 0x1.3731c0819ede4p-38, -0x1.f0caba642e0ebp-29, 0x1.6166f9d7a1c7fp-19, // subdomain 0:11 (0.16796875..0.171875)
                -0x1.2994e40bdcc8dp-87, 0x1.5208043f8ae32p-78, -0x1.6c285e8c598ccp-68, 0x1.71ae5f2e0ec43p-58, -0x1.5ec4df1e653b7p-48, 0x1.33bf63c6f12b5p-38, -0x1.ebf4dbbe1102ap-29, 0x1.5e816be03625cp-19, // subdomain 0:12 (0.171875..0.17578125)
                -0x1.24beeef2d32adp-87, 0x1.4cd7cf8517a0fp-78, -0x1.66eaa77e0c3a2p-68, 0x1.6cbd035fc6152p-58, -0x1.5a7742a08fec2p-48, 0x1.3057dbb39b1cfp-38, -0x1.e72cb0d51d8f4p-29, 0x1.5ba3146a8b7d7p-19, // subdomain 0:13 (0.17578125..0.1796875)
                -0x1.1ffe24a733d4ap-87, 0x1.47bd2dd52857cp-78, -0x1.61c185f6f068ap-68, 0x1.67dddb75b3775p-58, -0x1.56385ed64c52ap-48, 0x1.2cfb03353ec23p-38, -0x1.e2720ea108ee1p-29, 0x1.58cbdf095ed1dp-19, // subdomain 0:14 (0.1796875..0.18359375)
                -0x1.1b5224209bb99p-87, 0x1.42b7c0ca20522p-78, -0x1.5caca46853cadp-68, 0x1.63109ffa9b435p-58, -0x1.5207fd8f68832p-48, 0x1.29a8b5c204165p-38, -0x1.ddc4caacc0107p-29, 0x1.55fbb78f8ae7dp-19, // subdomain 0:15 (0.18359375..0.1875)
                -0x1.16ba8e26aba9cp-87, 0x1.3dc72baef1342p-78, -0x1.57abaeb9b6307p-68, 0x1.5e550aa224a0cp-58, -0x1.4de5e97053a90p-48, 0x1.2660cf56810eap-38, -0x1.d924bb1449ba9p-29, 0x1.53328a0f2d0a8p-19, // subdomain 0:16 (0.1875..0.19140625)
                -0x1.123705480535ep-87, 0x1.38eb137707426p-78, -0x1.52be524217845p-68, 0x1.59aad643c796dp-58, -0x1.49d1edeeaaff5p-48, 0x1.23232c73aa5a6p-38, -0x1.d491b682b0ed3p-29, 0x1.507042d8ccf55p-19, // subdomain 0:17 (0.19140625..0.1953125)
                -0x1.0dc72dd176076p-87, 0x1.34231eb65dce2p-78, -0x1.4de43dc165c8bp-68, 0x1.5511bed5d2540p-58, -0x1.45cbd74dd5dc0p-48, 0x1.1fefaa1cccb00p-38, -0x1.d00b942ff77abp-29, 0x1.4db4ce7a87d8ep-19, // subdomain 0:18 (0.1953125..0.19921875)
                -0x1.096aadc550b74p-87, 0x1.2f6ef599baf77p-78, -0x1.491d215a0a731p-68, 0x1.50898168853d9p-58, -0x1.41d3729bb0890p-48, 0x1.1cc625d58e812p-38, -0x1.cb922bdf10a66p-29, 0x1.4b0019bf3e6dep-19, // subdomain 0:19 (0.19921875..0.203125)
                -0x1.05212cd2f22edp-87, 0x1.2ace41df11fb2p-78, -0x1.4468ae8a968b5p-68, 0x1.4c11dc21455a3p-58, -0x1.3de88dad45a77p-48, 0x1.19a67d9ff9f18p-38, -0x1.c72555dbe3ac8p-29, 0x1.485211adc609fp-19, // subdomain 0:20 (0.203125..0.20703125)
                -0x1.00ea544e72a71p-87, 0x1.2640aece0b462p-78, -0x1.3fc698278cf5ep-68, 0x1.47aa8e35e4ad6p-58, -0x1.3a0af71b95dabp-48, 0x1.16908ffa8eec8p-38, -0x1.c2c4eaf95611ep-29, 0x1.45aaa3881caaep-19, // subdomain 0:21 (0.20703125..0.2109375)
                -0x1.f98b9e5102bc2p-88, 0x1.21c5e930b18bcp-78, -0x1.3b3692554a534p-68, 0x1.435357e8001b1p-58, -0x1.363a7e406d79fp-48, 0x1.13843bde5d360p-38, -0x1.be70c48f5d9b4p-29, 0x1.4309bcca9fecdp-19, // subdomain 0:22 (0.2109375..0.21484375)
                -0x1.f16693ccd43d1p-88, 0x1.1d5d9f4c4320ep-78, -0x1.36b8528209e27p-68, 0x1.3f0bfa80726c8p-58, -0x1.3276f33348083p-48, 0x1.108160bd26537p-38, -0x1.ba28bc7919ccbp-29, 0x1.406f4b2b46df6p-19, // subdomain 0:23 (0.21484375..0.21875)
                -0x1.e964e5348778ep-88, 0x1.190780da26e1fp-78, -0x1.324b8f6006d69p-68, 0x1.3ad4384adc120p-58, -0x1.2ec026c641341p-48, 0x1.0d87de7f872b9p-38, -0x1.b5ecad12f4d20p-29, 0x1.3ddb3c98deafcp-19, // subdomain 0:24 (0.21875..0.22265625)
                -0x1.e185f1b68d593p-88, 0x1.14c33f0103ed0p-78, -0x1.2df000dfb992ep-68, 0x1.36abd4913f3fbp-58, -0x1.2b15ea831321ep-48, 0x1.0a979583293aep-38, -0x1.b1bc7138cbb13p-29, 0x1.3b4d7f3a4a1b0p-19, // subdomain 0:25 (0.22265625..0.2265625)
                -0x1.d9c91b7c95c50p-88, 0x1.10908c4dfb7dap-78, -0x1.29a5602a30434p-68, 0x1.32929397b002fp-58, -0x1.277810a821c2fp-48, 0x1.07b06698fb3bfp-38, -0x1.ad97e4441da93p-29, 0x1.38c6016dc39f8p-19, // subdomain 0:26 (0.2265625..0.23046875)
                -0x1.d22dc79cd8486p-88, 0x1.0c6f1cae043b3p-78, -0x1.256b679b824a7p-68, 0x1.2e883a9817f3ap-58, -0x1.23e66c2592ff8p-48, 0x1.04d2330371236p-38, -0x1.a97ee20a429f7p-29, 0x1.3644b1c822617p-19, // subdomain 0:27 (0.23046875..0.234375)
                -0x1.cab35e0ba810ap-88, 0x1.085ea567664f6p-78, -0x1.2141d2bd5e00dp-68, 0x1.2a8c8fbe0d227p-58, -0x1.2060d09a73787p-48, 0x1.01fcdc74cb4fbp-38, -0x1.a57146daa87f1p-29, 0x1.33c97f1421b8bp-19, // subdomain 0:28 (0.234375..0.23828125)
                -0x1.c359498d51ae3p-88, 0x1.045edd13579bap-78, -0x1.1d285e41b0429p-68, 0x1.269f5a22bbea0p-58, -0x1.1ce71251e7a7bp-48, 0x1.fe608a1ac99a9p-39, -0x1.a16eef7d176cfp-29, 0x1.31545851ab5dap-19, // subdomain 0:29 (0.23828125..0.2421875)
                -0x1.bc1ef7a841156p-88, 0x1.006f7b97b76aap-78, -0x1.191ec7fd654dcp-68, 0x1.22c061c8e3496p-58, -0x1.197906406927cp-48, 0x1.f8d89eb411226p-39, -0x1.9d77b92ffcb4bp-29, 0x1.2ee52cb5242a4p-19, // subdomain 0:30 (0.2421875..0.24609375)
                -0x1.b503d8976e604p-88, 0x1.f9207441d1f85p-79, -0x1.1524cee34273bp-68, 0x1.1eef6f98e36fdp-58, -0x1.161682010fe9fp-48, 0x1.f361bca49b15ep-39, -0x1.998b81a6bc538p-29, 0x1.2c7beba6bb665p-19, // subdomain 0:31 (0.24609375..0.25)
                -0x1.aa95048b51e3bp-79, 0x1.edbead8ca0708p-71, -0x1.0f4a3a5464eacp-61, 0x1.194fd986332dep-52, -0x1.1117ff223c935p-43, 0x1.eb4ede1124bd6p-35, -0x1.93bd821bdcc39p-26, 0x1.28e8fe0d1be58p-17, // subdomain 1:0 (0.25..0.2578125)
                -0x1.9d1348576514dp-79, 0x1.defc8eac46f90p-71, -0x1.07b0f3e4543ccp-61, 0x1.11ff899944dd8p-52, -0x1.0a962b91443d9p-43, 0x1.e0c49bcd11481p-35, -0x1.8c2570d8cfd53p-26, 0x1.2439342ae7d67p-17, // subdomain 1:1 (0.2578125..0.265625)
                -0x1.9004bea5f7f3bp-79, 0x1.d0b1f60180194p-71, -0x1.0051c549aff83p-61, 0x1.0ae39e49404cap-52, -0x1.043f9be6afe1dp-43, 0x1.d67a927331918p-35, -0x1.84b7079cc7aa1p-26, 0x1.1f9ff3c039b13p-17, // subdomain 1:2 (0.265625..0.2734375)
                -0x1.836562917bad1p-79, 0x1.c2dae8f7e4f4ep-71, -0x1.f255b0de2b259p-62, 0x1.03fa87650422ep-52, -0x1.fc2634fba55f5p-44, 0x1.cc6f17708dfebp-35, -0x1.7d7148dd49831p-26, 0x1.1b1cc152f6b10p-17, // subdomain 1:3 (0.2734375..0.28125)
                -0x1.773154a618738p-79, 0x1.b573907d9022ep-71, -0x1.e474cdc3d3e43p-62, 0x1.fa8582bdffd42p-53, -0x1.f01ef5da0a64cp-44, 0x1.c2a08c1c37683p-35, -0x1.76533da233a9fp-26, 0x1.16af2457bc11dp-17, // subdomain 1:4 (0.28125..0.2890625)
                -0x1.6b64d97702037p-79, 0x1.a87837b888f7bp-71, -0x1.d6fd728dd3323p-62, 0x1.ed75a9c356b87p-53, -0x1.e467337628dc3p-44, 0x1.b90d5d5c6fd35p-35, -0x1.6f5bf556cbf6cp-26, 0x1.1256a71e6eca4p-17, // subdomain 1:5 (0.2890625..0.296875)
                -0x1.5ffc58420596cp-79, 0x1.9be54ac8ae64fp-71, -0x1.c9ec4da793ed2p-62, 0x1.e0c2accf77d54p-53, -0x1.d8fcb82e55c81p-44, 0x1.afb4034eb05e7p-35, -0x1.688a859c33ef7p-26, 0x1.0e12d6bf55febp-17, // subdomain 1:6 (0.296875..0.3046875)
                -0x1.54f459a0bd970p-79, 0x1.8fb75595a406bp-71, -0x1.bd3e299b75a64p-62, 0x1.d469cbafeb629p-53, -0x1.cddd5f263ae5cp-44, 0x1.a69300f2740f5p-35, -0x1.61de0a1d373e5p-26, 0x1.09e34308bc0dcp-17, // subdomain 1:7 (0.3046875..0.3125)
                -0x1.4a498646dd8b7p-79, 0x1.83eb02a8490d0p-71, -0x1.b0efec1931916p-62, 0x1.c8685c48ecb9bp-53, -0x1.c30713c00da56p-44, 0x1.9da8e3d6b0054p-35, -0x1.5b55a463699a9p-26, 0x1.05c77e6d12379p-17, // subdomain 1:8 (0.3125..0.3203125)
                -0x1.3ff8a5cd0eaafp-79, 0x1.787d1a0f458c7p-71, -0x1.a4fe9505470a2p-62, 0x1.bcbbc9dab3945p-53, -0x1.b877d11a3a656p-44, 0x1.94f443c9e2655p-35, -0x1.54f07bad899b4p-26, 0x1.01bf1df193015p-17, // subdomain 1:9 (0.3203125..0.328125)
                -0x1.35fe9d87db216p-79, 0x1.6d6a804e3f6e8p-71, -0x1.99673d9118b98p-62, 0x1.b161944d32bd1p-53, -0x1.ae2da1915c780p-44, 0x1.8c73c28ca21e3p-35, -0x1.4eadbcc71e45dp-26, 0x1.fb93723abf2b5p-18, // subdomain 1:10 (0.328125..0.3359375)
                -0x1.2c586f6a2b840p-79, 0x1.62b035574cb00p-71, -0x1.8e27175b66cd4p-62, 0x1.a6574f8210834p-53, -0x1.a4269e465ce57p-44, 0x1.84260b869a63cp-35, -0x1.488c99e1457f0p-26, 0x1.f3cdd3d22901ep-18, // subdomain 1:11 (0.3359375..0.34375)
                -0x1.230338f2df1f5p-79, 0x1.584b538e3ce7ep-71, -0x1.833b6b98d5f0bp-62, 0x1.9b9aa2ac9f802p-53, -0x1.9a60eea8933f1p-44, 0x1.7c09d37ddd7c6p-35, -0x1.428c4a6ca9d8bp-26, 0x1.ec2c995da488bp-18, // subdomain 1:12 (0.34375..0.3515625)
                -0x1.19fc32250cbe5p-79, 0x1.4e390ed5582aap-71, -0x1.78a19a4435d81p-62, 0x1.912947afa0659p-53, -0x1.90dac803c526cp-44, 0x1.741dd8507b2f4p-35, -0x1.3cac0af496869p-26, 0x1.e4af00356cb11p-18, // subdomain 1:13 (0.3515625..0.359375)
                -0x1.1140ac8a7e4bap-79, 0x1.4476b3a3354e2p-71, -0x1.6e5719563d2f1p-62, 0x1.87010a80979d1p-53, -0x1.87926d11e276dp-44, 0x1.6c60e0b047ceep-35, -0x1.36eb1cfb20923p-26, 0x1.dd544a3782a40p-18, // subdomain 1:14 (0.359375..0.3671875)
                -0x1.08ce123ffe38ep-79, 0x1.3b01a6214d5a0p-71, -0x1.645974047399dp-62, 0x1.7d1fc890837dap-53, -0x1.7e862d905d2ccp-44, 0x1.64d1bbe0c184fp-35, -0x1.3148c6d65caf8p-26, 0x1.d61bbdaa4dea1p-18, // subdomain 1:15 (0.3671875..0.375)
                -0x1.00a1e50b11f21p-79, 0x1.31d76152f58abp-71, -0x1.5aa64a07053ddp-62, 0x1.73837039c1c60p-53, -0x1.75b465d8fd689p-44, 0x1.5d6f417702227p-35, -0x1.2bc4538e995e3p-26, 0x1.cf04a5200a3f1p-18, // subdomain 1:16 (0.375..0.3828125)
                -0x1.f1737af1819aep-80, 0x1.28f576446cc9cp-71, -0x1.513b4ee53dec1p-62, 0x1.6a2a0032f4e20p-53, -0x1.6d1b7e7e12f54p-44, 0x1.5638511bac562p-35, -0x1.265d12bd9546bp-26, 0x1.c80e4f5afd107p-18, // subdomain 1:17 (0.3828125..0.390625)
                -0x1.e226940a15ccbp-80, 0x1.20598b41accc4p-71, -0x1.481649486c938p-62, 0x1.61118706bb433p-53, -0x1.64b9ebe9e6e44p-44, 0x1.4f2bd24ec3beep-35, -0x1.2112586eaa07ap-26, 0x1.c1380f326cd3ep-18, // subdomain 1:18 (0.390625..0.3984375)
                -0x1.d3589c956e443p-80, 0x1.18015b14b23bbp-71, -0x1.3f351254f2129p-62, 0x1.583822900ccaep-53, -0x1.5c8e2e0150d58p-44, 0x1.4848b42d5fe13p-35, -0x1.1be37cffe3e80p-26, 0x1.ba813b7854938p-18, // subdomain 1:19 (0.3984375..0.40625)
                -0x1.c5054476468c5p-80, 0x1.0feab44af25a0p-71, -0x1.3695950940d8dp-62, 0x1.4f9bff7b13e25p-53, -0x1.5496cfc954871p-44, 0x1.418ded392a923p-35, -0x1.16cfdd03ff210p-26, 0x1.b3e92edfce36cp-18, // subdomain 1:20 (0.40625..0.4140625)
                -0x1.b72862bb890ffp-80, 0x1.08137881b77f9p-71, -0x1.2e35cda295096p-62, 0x1.473b58ca59733p-53, -0x1.4cd2670fad417p-44, 0x1.3afa7b219af1dp-35, -0x1.11d6d92543b10p-26, 0x1.ad6f47e42e3cdp-18, // subdomain 1:21 (0.4140625..0.421875)
                -0x1.a9bdf42e14b0bp-80, 0x1.00799bb920a18p-71, -0x1.2613c9072ef8dp-62, 0x1.3f1477601c728p-53, -0x1.453f94162da63p-44, 0x1.348d628ede8e0p-35, -0x1.0cf7d60938e0ep-26, 0x1.a712e8b0cbc5bp-18, // subdomain 1:22 (0.421875..0.4296875)
                -0x1.9cc219ecb0c59p-80, 0x1.f236475b05811p-72, -0x1.1e2da435dbf45p-62, 0x1.3725b18baf251p-53, -0x1.3ddd0140db4d5p-44, 0x1.2e45aeee62c7ap-35, -0x1.08323c3529dffp-26, 0x1.a0d377096fefbp-18, // subdomain 1:23 (0.4296875..0.4375)
                -0x1.903118159f986p-80, 0x1.e3ec4e6ddb6a0p-72, -0x1.16818bba9b407p-62, 0x1.2f6d6a9ab5aa3p-53, -0x1.36a962c6ae7dcp-44, 0x1.28227240f1094p-35, -0x1.038577f375108p-26, 0x1.9ab05c3367b58p-18, // subdomain 1:24 (0.4375..0.4453125)
                -0x1.8407547d455c6p-80, 0x1.d6119b5d30902p-72, -0x1.0f0dbb282f249p-62, 0x1.27ea126e22af7p-53, -0x1.2fa37664df228p-44, 0x1.2222c4ea50d4ep-35, -0x1.fde1f2733fa6cp-27, 0x1.94a904df33a0dp-18, // subdomain 1:25 (0.4453125..0.453125)
                -0x1.784155715f3e5p-80, 0x1.c8a29cb3f3d26p-72, -0x1.07d07c966bb3bp-62, 0x1.209a2512d07e2p-53, -0x1.28ca0314a8ee1p-44, 0x1.1c45c582631dcp-35, -0x1.f4e8671e51a8dp-27, 0x1.8ebce112d0c6dp-18, // subdomain 1:26 (0.453125..0.4609375)
                -0x1.6cdbc0884b93ap-80, 0x1.bb9bdfb6b52f7p-72, -0x1.00c8282516c2dp-62, 0x1.197c2a5d95ccfp-53, -0x1.221bd8c371644p-44, 0x1.168a98a7aad05p-35, -0x1.ec1d3bea26171p-27, 0x1.88eb641496aabp-18, // subdomain 1:27 (0.4609375..0.46875)
                -0x1.61d3597be95dcp-80, 0x1.aefa0f4f478e5p-72, -0x1.f3e647065c667p-63, 0x1.128eb58ab6f63p-53, -0x1.1b97d00d3b429p-44, 0x1.10f068d336d1fp-35, -0x1.e37f659c47a9ap-27, 0x1.83340456a5c18p-18, // subdomain 1:28 (0.46875..0.4765625)
                -0x1.5725010f9a4a8p-80, 0x1.a2b9f3027e3b8p-72, -0x1.e69fc2f4d6ae5p-63, 0x1.0bd064e09556cp-53, -0x1.153cc9f953780p-44, 0x1.0b76662de227bp-35, -0x1.db0ddf8f9d608p-27, 0x1.7d963b62e2765p-18, // subdomain 1:29 (0.4765625..0.484375)
                -0x1.4ccdb400f71c3p-80, 0x1.96d86defa3e04p-72, -0x1.d9b9c2fdb78d4p-63, 0x1.053fe1557fb52p-53, -0x1.0f09afb9248f2p-44, 0x1.061bc666df541p-35, -0x1.d2c7ab8792ceep-27, 0x1.781185c778abap-18, // subdomain 1:30 (0.484375..0.4921875)
                -0x1.42ca8a02caf09p-80, 0x1.8b527dd94c95ep-72, -0x1.cd315e84d2fc5p-63, 0x1.fdb7bc710f44ap-54, -0x1.08fd72691e190p-44, 0x1.00dfc48b845adp-35, -0x1.caabd18487013p-27, 0x1.72a56303e3d7ap-18, // subdomain 1:31 (0.4921875..0.5)
                -0x1.345ee3543075bp-71, 0x1.7ab0e6db222b1p-64, -0x1.bb0c202f399abp-56, 0x1.eb2d1654a35cfp-48, -0x1.0031b6f9f969bp-39, 0x1.f27b0fabcb6f6p-32, -0x1.becf63b1d12e0p-24, 0x1.6ab030260bd40p-16, // subdomain 2:0 (0.5..0.515625)
                -0x1.222fc2964d4bdp-71, 0x1.65a79bfc1307ap-64, -0x1.a40b61677b864p-56, 0x1.d394d59a38787p-48, -0x1.e9ec24740392ap-40, 0x1.dee86dbf4d613p-32, -0x1.af8a6b80616c7p-24, 0x1.60657102402fap-16, // subdomain 2:1 (0.515625..0.53125)
                -0x1.1121c0ee31e79p-71, 0x1.51db5a93ae4b3p-64, -0x1.8e512ee0f911cp-56, 0x1.bd359c868c6ccp-48, -0x1.d48895f369533p-40, 0x1.cc3105d0b2901p-32, -0x1.a0de944841e21p-24, 0x1.567480dfd309ap-16, // subdomain 2:2 (0.53125..0.546875)
                -0x1.01220ca1dc4f2p-71, 0x1.3f3879c69f470p-64, -0x1.79ca576281678p-56, 0x1.a7fe171323467p-48, -0x1.c02a7f1b284dap-40, 0x1.ba4a5a8ffa0a4p-32, -0x1.92c52e8dd8566p-24, 0x1.4cd9dd315404dp-16, // subdomain 2:3 (0.546875..0.5625)
                -0x1.e43e36889c769p-72, 0x1.2dac97d383338p-64, -0x1.6664d9e36a641p-56, 0x1.93ddf583b1990p-48, -0x1.acc4665a32f12p-40, 0x1.a92a795380557p-32, -0x1.8537dc9030113p-24, 0x1.43922a8edddf5p-16, // subdomain 2:4 (0.5625..0.578125)
                -0x1.c8112529e4098p-72, 0x1.1d2683e23dd27p-64, -0x1.540fd1ce08bf0p-56, 0x1.80c5dc4f1b26fp-48, -0x1.9a4990138dff3p-40, 0x1.98c7f271fe4ebp-32, -0x1.78308e12adbe8p-24, 0x1.3a9a32d8780dep-16, // subdomain 2:5 (0.578125..0.59375)
                -0x1.ad9e6621adc19p-72, 0x1.0d96296425e5cp-64, -0x1.42bb64987bf51p-56, 0x1.6ea7551358325p-48, -0x1.88adf38734974p-40, 0x1.8919d20f0a756p-32, -0x1.6ba97c6224746p-24, 0x1.31eee3710a317p-16, // subdomain 2:6 (0.59375..0.609375)
                -0x1.94c98829c4e83p-72, 0x1.fdd8f9cf687bdp-65, -0x1.3258b09baec7ep-56, 0x1.5d74c08233aa7p-48, -0x1.77e63069a7068p-40, 0x1.7a1799530e3f4p-32, -0x1.5f9d2691d4527p-24, 0x1.298d4b9087a4ep-16, // subdomain 2:7 (0.609375..0.625)
                -0x1.7d7802e5507bfp-72, 0x1.e236d489311e2p-65, -0x1.22d9bd12306dfp-56, 0x1.4d21493619a57p-48, -0x1.67e7851ee01e5p-40, 0x1.6bb9380815c75p-32, -0x1.54064ded11cbep-24, 0x1.21729aac0d84fp-16, // subdomain 2:8 (0.625..0.640625)
                -0x1.679114b7e99aap-72, 0x1.c82b880e467a0p-65, -0x1.14316b2a2bb9ap-56, 0x1.3da0d75f5b0e2p-48, -0x1.58a7c58e2f120p-40, 0x1.5df706954cfc3p-32, -0x1.48dff29a87c2bp-24, 0x1.199c1ef2c0821p-16, // subdomain 2:9 (0.640625..0.65625)
                -0x1.52fda3168ee96p-72, 0x1.af9e6807acfabp-65, -0x1.06536817432f8p-56, 0x1.2ee8053b62fa7p-48, -0x1.4a1d528528809p-40, 0x1.50c9c05362f52p-32, -0x1.3e25506e36674p-24, 0x1.120743de498abp-16, // subdomain 2:10 (0.65625..0.671875)
                -0x1.3fa81d22eaab3p-72, 0x1.987859fb5d9b0p-65, -0x1.f2684024fb7a1p-57, 0x1.20ec14485d5e1p-48, -0x1.3c3f11a089517p-40, 0x1.442a7e346ac6bp-32, -0x1.33d1dbe77d8a1p-24, 0x1.0ab190d5e10d2p-16, // subdomain 2:11 (0.671875..0.6875)
                -0x1.2d7c60650f8d0p-72, 0x1.82a3ba969b997p-65, -0x1.d991646f73611p-57, 0x1.13a2e328bf485p-48, -0x1.2f0465b28259fp-40, 0x1.3812b1ba2990ap-32, -0x1.29e13f58aaee4p-24, 0x1.0398a7e2e8350p-16, // subdomain 2:12 (0.6875..0.703125)
                -0x1.1c679f7b3359ap-72, 0x1.6e0c44d15ba43p-65, -0x1.c20dca6299b26p-57, 0x1.0702e42afffb9p-48, -0x1.2265279e7534ep-40, 0x1.2c7c203613f75p-32, -0x1.204f5835ac5ebp-24, 0x1.f97488ec3c826p-17, // subdomain 2:13 (0.703125..0.71875)
                -0x1.0c584a9a0cf04p-72, 0x1.5a9efac8c53ddp-65, -0x1.abca394726687p-57, 0x1.f60628d542780p-49, -0x1.16599fa2b212fp-40, 0x1.2160de4e8a3aep-32, -0x1.17183487aa001p-24, 0x1.ec2874791bbebp-17, // subdomain 2:14 (0.71875..0.734375)
                -0x1.fa7bf376b346fp-73, 0x1.484a10315d47fp-65, -0x1.96b49e9666076p-57, 0x1.df35e700cac77p-49, -0x1.0ada7f09444e6p-40, 0x1.16bb4bc52a031p-32, -0x1.0e3810836fbbep-24, 0x1.df48e806b6f04p-17, // subdomain 2:15 (0.734375..0.75)
                -0x1.de12b0b565865p-73, 0x1.36fcd643bd3e8p-65, -0x1.82bbfb839e78fp-57, 0x1.c984f770886efp-49, -0x1.ffc1b472a1c81p-41, 0x1.0c860f7a501a0p-32, -0x1.05ab5440acd9fp-24, 0x1.d2d1e5535316fp-17, // subdomain 2:16 (0.75..0.765625)
                -0x1.c358254696775p-73, 0x1.26a7a908e19aep-65, -0x1.6fd053bb82bdap-57, 0x1.b4e43588aa511p-49, -0x1.eacc4645ee6dep-41, 0x1.02bc13aa2309dp-32, -0x1.fadd23205e19cp-25, 0x1.c6bf983fa76dap-17, // subdomain 2:17 (0.765625..0.78125)
                -0x1.aa3187e6d6893p-73, 0x1.173bdded16b3dp-65, -0x1.5de29d4486c97p-57, 0x1.a14555c8ff240p-49, -0x1.d6c84800149e9p-41, 0x1.f2b104bf74501p-33, -0x1.eafd03deb505ap-25, 0x1.bb0e54e697931p-17, // subdomain 2:18 (0.78125..0.796875)
                -0x1.9285cd3922113p-73, 0x1.08abb382618b4p-65, -0x1.4ce4b16c65837p-57, 0x1.8e9ad8f6a0ea0p-49, -0x1.c3a9f4c8af327p-41, 0x1.e0ad843656bc1p-33, -0x1.dbb0092e761fcp-25, 0x1.afba95ccdce4ap-17, // subdomain 2:19 (0.796875..0.8125)
                -0x1.7c3d899209490p-73, 0x1.f5d484b9e8c3ap-66, -0x1.3cc93eb08f3e3p-57, 0x1.7cd80011b4185p-49, -0x1.b1662531110e8p-41, 0x1.cf64e54efa564p-33, -0x1.ccf0394b46d12p-25, 0x1.a4c0fa275d36cp-17, // subdomain 2:20 (0.8125..0.828125)
                -0x1.6742d4e2de679p-73, 0x1.dbd6ddeba0c52p-66, -0x1.2d83bb9081448p-57, 0x1.6bf0c109fa01dp-49, -0x1.9ff24675fdb62p-41, 0x1.beced43b546e2p-33, -0x1.beb7debdcd66dp-25, 0x1.9a1e4436ff10ep-17, // subdomain 2:21 (0.828125..0.84375)
                -0x1.5381309d45ca1p-73, 0x1.c347b90251b80p-66, -0x1.1f085a3a39381p-57, 0x1.5bd9bc25db9a0p-49, -0x1.8f445244b03f4p-41, 0x1.aee364242112bp-33, -0x1.b101850d9d32ap-25, 0x1.8fcf57b8dc77cp-17, // subdomain 2:22 (0.84375..0.859375)
                -0x1.40e56f6a8035bp-73, 0x1.ac11c14ee74aep-66, -0x1.114bfd0212d6ap-57, 0x1.4c8832106706bp-49, -0x1.7f52c6fb18f61p-41, 0x1.9f9b09d3a0a50p-33, -0x1.a3c7f59f0907ap-25, 0x1.85d13869c53cbp-17, // subdomain 2:23 (0.859375..0.875)
                -0x1.2f5d9e9487f87p-73, 0x1.9620f0229f51bp-66, -0x1.04442b9865ed0p-57, 0x1.3df1fa8380e16p-49, -0x1.7014a05bbf864p-41, 0x1.90ee96ab8ab1ap-33, -0x1.970634ba6ff14p-25, 0x1.7c21089c10c23p-17, // subdomain 2:24 (0.875..0.890625)
                -0x1.1ed8f101a3d9ap-73, 0x1.8162776f202dep-66, -0x1.efce11e05cbf3p-58, 0x1.300d7b8442996p-49, -0x1.618150ae3418bp-41, 0x1.82d733e1d5742p-33, -0x1.8ab77ebabc2d9p-25, 0x1.72bc07decd427p-17, // subdomain 2:25 (0.890625..0.90625)
                -0x1.0f47aba56abaap-73, 0x1.6dc4add394c18p-66, -0x1.d8569395d552bp-58, 0x1.22d1a1282b194p-49, -0x1.5390ba456e0dfp-41, 0x1.754e5e0035f33p-33, -0x1.7ed74560ee617p-25, 0x1.699f91b567ec1p-17, // subdomain 2:26 (0.90625..0.921875)
                -0x1.009b134c5c7abp-73, 0x1.5b36fbfd80c88p-66, -0x1.c21057bd5eb09p-58, 0x1.1635d5da66959p-49, -0x1.463b2965e3771p-41, 0x1.684de0a27ccd1p-33, -0x1.73612d4aaedb6p-25, 0x1.60c91c5ef5ae7p-17, // subdomain 2:27 (0.921875..0.9375)
                -0x1.e58ab75267a7cp-74, 0x1.49a9cb45d1ac3p-66, -0x1.aceadb1862298p-58, 0x1.0a31fb19071abp-49, -0x1.39794e8589ad1p-41, 0x1.5bcfd2702fdc1p-33, -0x1.68510b89fca49p-25, 0x1.583637ac5057cp-17, // subdomain 2:28 (0.9375..0.953125)
                -0x1.cb732f19d9adep-74, 0x1.390e75745e932p-66, -0x1.98d68c3677993p-58, 0x1.fd7cc54523b07p-50, -0x1.2d4438e04deb5p-41, 0x1.4fce914df91a6p-33, -0x1.5da2e35c2e3d5p-25, 0x1.4fe48be446c49p-17, // subdomain 2:29 (0.953125..0.96875)
                -0x1.b2d754685f02dp-74, 0x1.295735959938cp-66, -0x1.85c4bce96b61cp-58, 0x1.e7a790198c0ccp-50, -0x1.2195515bed130p-41, 0x1.4444bec3b7d72p-33, -0x1.5352e3fe9247bp-25, 0x1.47d1d8b529588p-17, // subdomain 2:30 (0.96875..0.984375)
                -0x1.9ba0732cc0f8bp-74, 0x1.1a7719cfa26a4p-66, -0x1.73a794a3d19b6p-58, 0x1.d2d6957bf487cp-50, -0x1.166655b6642bep-41, 0x1.392d3c9432a95p-33, -0x1.495d669f17342p-25, 0x1.3ffbf43305d77p-17, // subdomain 2:31 (0.984375..1.0)
                -0x1.7b43026c2f569p-65, 0x1.05a429e3a17fap-58, -0x1.5a2532f0b696dp-51, 0x1.b568c166e43b1p-44, -0x1.0682ca268c417p-36, 0x1.2955b33169bc3p-29, -0x1.3b0f4350ef3bep-22, 0x1.34a899c7ce572p-15, // subdomain 3:0 (1.0..1.03125)
                -0x1.543920f95dc27p-65, 0x1.d8bee17dd874ap-59, -0x1.3b0da62b5a487p-51, 0x1.9146a69b9e330p-44, -0x1.e5bc28643a2d2p-37, 0x1.159abd0f5105fp-29, -0x1.2919cb8fbb2d1p-22, 0x1.2650df14e479bp-15, // subdomain 3:1 (1.03125..1.0625)
                -0x1.316cd8ba21294p-65, 0x1.ab64f95962bb6p-59, -0x1.1ef442fb3702ap-51, 0x1.70608992b5dc2p-44, -0x1.c1ae77855f68dp-37, 0x1.035839d03e4f6p-29, -0x1.1854147ee6cdbp-22, 0x1.18c96ff487903p-15, // subdomain 3:2 (1.0625..1.09375)
                -0x1.126261f01473fp-65, 0x1.82aad44d37b1ap-59, -0x1.058aaa77f175fp-51, 0x1.5266f9ee6551dp-44, -0x1.a092ac4ee0cdbp-37, 0x1.e4de970a56aadp-30, -0x1.08a791b7fea62p-22, 0x1.0c049896a4e8fp-15, // subdomain 3:3 (1.09375..1.125)
                -0x1.ed5969009c2ccp-66, 0x1.5e110423056c5p-59, -0x1.dd15e8388eabdp-52, 0x1.3712a110261d0p-44, -0x1.822600d916cf9p-37, 0x1.c587a54b0d727p-30, -0x1.f3ff1cc45a76ep-23, 0x1.ffeb510d824ddp-16, // subdomain 3:4 (1.125..1.15625)
                -0x1.bbd7672bcb01fp-66, 0x1.3d268235470bbp-59, -0x1.b36d7055d2abbp-52, 0x1.1e2362d604172p-44, -0x1.662c06f9439f8p-37, 0x1.a878563923162p-30, -0x1.d89206d9f43cap-23, 0x1.e921bb2046a3fp-16, // subdomain 3:5 (1.15625..1.1875)
                -0x1.8f952b20ccd05p-66, 0x1.1f86fbb009be4p-59, -0x1.8dac6e5c2ad0dp-52, 0x1.075f979551cd6p-44, -0x1.4c6e043d7f252p-37, 0x1.8d81e03f13fedp-30, -0x1.bee4e385cba31p-23, 0x1.d396a0fe067d6p-16, // subdomain 3:6 (1.1875..1.21875)
                -0x1.67fd006b09b7dp-66, 0x1.04d9537174977p-59, -0x1.6b6f35f5134f8p-52, 0x1.e526b8ae48462p-45, -0x1.34ba5f68783bap-37, 0x1.74798c2585795p-30, -0x1.a6d7801a84b79p-23, 0x1.bf35c72320a44p-16, // subdomain 3:7 (1.21875..1.25)
                -0x1.448aa470a9fcdp-66, 0x1.d99ca560b3007p-60, -0x1.4c5c9c28574eep-52, 0x1.bf1fed71ec2bdp-45, -0x1.1ee41d76f34e1p-37, 0x1.5d3853b156629p-30, -0x1.904c375b3779ap-23, 0x1.abec64e0644afp-16, // subdomain 3:8 (1.25..1.28125)
                -0x1.24c928363c577p-66, 0x1.ae3f04ee07f5dp-60, -0x1.3024cd71cc48ap-52, 0x1.9c56944b9d903p-45, -0x1.0ac26c73c751bp-37, 0x1.479a89f07b5a0p-30, -0x1.7b27b93579487p-23, 0x1.99a90710ea50fp-16, // subdomain 3:9 (1.28125..1.3125)
                -0x1.08511682e873ap-66, 0x1.871c519176f87p-60, -0x1.1680465dd6b47p-52, 0x1.7c7eb8c7a3b98p-45, -0x1.f0607525873e8p-38, 0x1.337f8c31210f3p-30, -0x1.6750d7afff0d5p-23, 0x1.885b75536fe5bp-16, // subdomain 3:10 (1.3125..1.34375)
                -0x1.dd8daaee2b23cp-67, 0x1.63c4c9eaf054cp-60, -0x1.fe5dd5024b339p-53, 0x1.5f53cf6147a86p-45, -0x1.ce17b07098fb6p-38, 0x1.20c97abb574a7p-30, -0x1.54b0588a95ba0p-23, 0x1.77f4998bdd692p-16, // subdomain 3:11 (1.34375..1.375)
                -0x1.afb277b7c5b2fp-67, 0x1.43d4ce69d42dap-60, -0x1.d3ee6c4c08546p-53, 0x1.4497f1bb61705p-45, -0x1.ae6d4776c7d78p-38, 0x1.0f5cf87efba38p-30, -0x1.4330cb19dabdfp-23, 0x1.6866697777e97p-16, // subdomain 3:12 (1.375..1.40625)
                -0x1.8680a6f6b8640p-67, 0x1.26f38046dffa5p-60, -0x1.ad4b0efc5363fp-53, 0x1.2c13301a62b7cp-45, -0x1.91297ba5e0efap-38, 0x1.fe41e1fc4431fp-31, -0x1.32be61f54ac23p-23, 0x1.59a3d22374325p-16, // subdomain 3:13 (1.40625..1.4375)
                -0x1.6178856c3f5f9p-67, 0x1.0cd18a577db05p-60, -0x1.8a16f91acee58p-53, 0x1.1592f5b54b579p-45, -0x1.76198b7bfbf9cp-38, 0x1.dffcc7a05b547p-31, -0x1.2346d0191643cp-23, 0x1.4ba0a51a51024p-16, // subdomain 3:14 (1.4375..1.46875)
                -0x1.4028b0a0d3927p-67, 0x1.ea501d24b272dp-61, -0x1.69fec6b3f9d75p-53, 0x1.00e97dbc4b370p-45, -0x1.5d0f38a3e29dfp-38, 0x1.c3c06a56f48f4p-31, -0x1.14b92926d8303p-23, 0x1.3e5187209197ap-16, // subdomain 3:15 (1.46875..1.5)
                -0x1.222c6a3564e9fp-67, 0x1.bf6f678b716ddp-61, -0x1.4cb77478e470fp-53, 0x1.dbdaae684213fp-46, -0x1.45e05a635907ap-38, 0x1.a96605e6518ecp-31, -0x1.0705c478ec987p-23, 0x1.31abe05d37af5p-16, // subdomain 3:16 (1.5..1.53125)
                -0x1.072a1f59f9bc2p-67, 0x1.988f9f24d2923p-61, -0x1.31fd7cd83a2e1p-53, 0x1.b8f1ebfb6011bp-46, -0x1.30667b0ba5db4p-38, 0x1.90c9fba912c65p-31, -0x1.f43c4587b5e01p-24, 0x1.25a5cdcdc0076p-16, // subdomain 3:17 (1.53125..1.5625)
                -0x1.dda43babea6b6p-68, 0x1.754b527a59ae1p-61, -0x1.19940e4065f17p-53, 0x1.98d49f259f1d4p-46, -0x1.1c7e7f4225093p-38, 0x1.79cb8c17fd2b0p-31, -0x1.dbe9ac1055508p-24, 0x1.1a3613e870f3cp-16, // subdomain 3:18 (1.5625..1.59375)
                -0x1.b1bae1d3d6156p-68, 0x1.5547a4ae21ee5p-61, -0x1.034457a13ef40p-53, 0x1.7b450459bb158p-46, -0x1.0a085617c32cbp-38, 0x1.644c96e85ef33p-31, -0x1.c4fad73a436c5p-24, 0x1.0f54125288ddap-16, // subdomain 3:19 (1.59375..1.625)
                -0x1.8a19c221f385ep-68, 0x1.3833253a63a83p-61, -0x1.ddb9d13ae88e7p-54, 0x1.600b0a5fc3048p-46, -0x1.f1cd6206fa8cdp-39, 0x1.503161184de58p-31, -0x1.af58b044ded46p-24, 0x1.04f7b89252dedp-16, // subdomain 3:20 (1.625..1.65625)
                -0x1.66502498807e9p-68, 0x1.1dc4c988282bap-61, -0x1.b8624654fa818p-54, 0x1.46f3c3ad73d80p-46, -0x1.d1fd85dbc15c6p-39, 0x1.3d6060535f95dp-31, -0x1.9aedca576a74ep-24, 0x1.f632f74ec2c94p-17, // subdomain 3:21 (1.65625..1.6875)
                -0x1.45f97a6354ba4p-68, 0x1.05bb045740f79p-61, -0x1.96317945442ebp-54, 0x1.2fd0e674e73f2p-46, -0x1.b4700b24e0c26p-39, 0x1.2bc20b29fa7f3p-31, -0x1.87a640fa95e5dp-24, 0x1.e36498ea60e60p-17, // subdomain 3:22 (1.6875..1.71875)
                -0x1.28bbfee610198p-68, 0x1.dfb5eee149033p-62, -0x1.76de9911f57d3p-54, 0x1.1a7859d215239p-46, -0x1.98f80920438aap-39, 0x1.1b40ada32defep-31, -0x1.756f996b24f8fp-24, 0x1.d1771ddc81088p-17, // subdomain 3:23 (1.71875..1.75)
                -0x1.0e4781fb81acbp-68, 0x1.b7df79164e42ap-62, -0x1.5a27cee842dc0p-54, 0x1.06c3ceadb70f0p-46, -0x1.7f6c5b2845c02p-39, 0x1.0bc841bd4f396p-31, -0x1.6438a6837b931p-24, 0x1.c05e23debab54p-17, // subdomain 3:24 (1.75..1.78125)
                -0x1.eca8acc2d5af9p-69, 0x1.939385568c007p-62, -0x1.3fd1898578462p-54, 0x1.e920c629837d0p-47, -0x1.67a74a3428aefp-39, 0x1.fa8c96f92f9c2p-32, -0x1.53f16f0035761p-24, 0x1.b00e0f787befdp-17, // subdomain 3:25 (1.78125..1.8125)
                -0x1.c144bfdbae116p-69, 0x1.727c7b367a891p-62, -0x1.27a5dbcb7de84p-54, 0x1.c77c9dcf17631p-47, -0x1.51863e961eb6cp-39, 0x1.df537061a0200p-32, -0x1.448b15e9ade83p-24, 0x1.a07bfe2000ff0p-17, // subdomain 3:26 (1.8125..1.84375)
                -0x1.99f07b706de54p-69, 0x1.544d4ac5062dep-62, -0x1.1173ec602d925p-54, 0x1.a86131c6d79c6p-47, -0x1.3ce97922d0c17p-39, 0x1.c5c58147ae838p-32, -0x1.35f7c4f182f7bp-24, 0x1.919db96a7363ep-17, // subdomain 3:27 (1.84375..1.875)
                -0x1.764519b63885fp-69, 0x1.38c087c8c14fbp-62, -0x1.fa1ee8eb9f60bp-55, 0x1.8b99a64b45633p-47, -0x1.29b3d304a4e89p-39, 0x1.adc59bd56010cp-32, -0x1.282a9897c32c2p-24, 0x1.8369ab34ce84ap-17, // subdomain 3:28 (1.875..1.90625)
                -0x1.55e681bc944c1p-69, 0x1.1f979e0be708dp-62, -0x1.d4a09824fd234p-55, 0x1.70f5bcbe5a685p-47, -0x1.17ca8390fa799p-39, 0x1.9738c133792a3p-32, -0x1.1b178defa0dd4p-24, 0x1.75d6d2b122dfep-17, // subdomain 3:29 (1.90625..1.9375)
                -0x1.38821e6bd6e9ep-69, 0x1.089a1bdadc0bbp-62, -0x1.b22404b9f003dp-55, 0x1.5849660e245c8p-47, -0x1.0714eb874be93p-39, 0x1.8205f4463a3c1p-32, -0x1.0eb371df4888ap-24, 0x1.68dcba459b652p-17, // subdomain 3:30 (1.9375..1.96875)
                -0x1.1dcdd729e5ab7p-69, 0x1.e72a205f6f3e0p-63, -0x1.92670140c2e4cp-55, 0x1.416c5fd85eefcp-47, -0x1.eef8ca63e71c4p-40, 0x1.6e16105786c9ap-32, -0x1.02f3d1b9da858p-24, 0x1.5c736e2c48839p-17, // subdomain 3:31 (1.96875..2.0)
                -0x1.f487fa6bf34dap-61, 0x1.aeca204d698a4p-55, -0x1.6763734f98344p-48, 0x1.2235e605eed8ap-41, -0x1.c401f3dfe74bcp-35, 0x1.525d3407e8af9p-28, -0x1.e4e744ecff94bp-22, 0x1.4ad4be02541f5p-15, // subdomain 4:0 (2.0..2.0625)
                -0x1.a409344393754p-61, 0x1.6e3c7d6cce811p-55, -0x1.35b6cd60cd4b8p-48, 0x1.fb5f0313e0369p-42, -0x1.911a2c2b57428p-35, 0x1.311583a189f1cp-28, -0x1.bcbab37ef56a9p-22, 0x1.350f1adcd7194p-15, // subdomain 4:1 (2.0625..2.125)
                -0x1.61402a0183d41p-61, 0x1.3801df2d99db7p-55, -0x1.0b70942b4c7ddp-48, 0x1.bc5b571d64bc2p-42, -0x1.648fed6515ad5p-35, 0x1.1386b63e574bcp-28, -0x1.987a369381f6cp-22, 0x1.21138efb1a061p-15, // subdomain 4:2 (2.125..2.1875)
                -0x1.29b794d957475p-61, 0x1.0a5a7a17bf9fep-55, -0x1.cec6f6c724a02p-49, 0x1.85e44119c7183p-42, -0x1.3d84b6ef1a81fp-35, 0x1.f274d33e0e14bp-29, -0x1.77b6464382936p-22, 0x1.0eb5b56134517p-15, // subdomain 4:3 (2.1875..2.25)
                -0x1.f6e0b8d5e084cp-62, 0x1.c7acca13a4355p-56, -0x1.9129fb0ecac85p-49, 0x1.56b96ac335517p-42, -0x1.1b3acc5699da5p-35, 0x1.c394990186698p-29, -0x1.5a0d1f0d5a658p-22, 0x1.fb9c1dab9cf2cp-16, // subdomain 4:4 (2.25..2.3125)
                -0x1.a994dfd0691cap-62, 0x1.868d5615a3308p-56, -0x1.5c68b17d6b6cfp-49, 0x1.2dcc975919640p-42, -0x1.fa201b9290b6ep-36, 0x1.99bb88d125e0fp-29, -0x1.3f28e66b6ad9cp-22, 0x1.dc72d03a2ae68p-16, // subdomain 4:5 (2.3125..2.375)
                -0x1.68e4ac2308054p-62, 0x1.4f628cd0e5a29p-56, -0x1.2f2671aad1acep-49, 0x1.0a3954935880cp-42, -0x1.c4f359357f64dp-36, 0x1.744fe34b7be3bp-29, -0x1.26be149fde732p-22, 0x1.bfb0a9a336cd7p-16, // subdomain 4:6 (2.375..2.4375)
                -0x1.32a4eac4ec85dp-62, 0x1.208edc19bdd68p-56, -0x1.083f760ddf256p-49, 0x1.d67c3eb4e06d4p-43, -0x1.96012a5d42116p-36, 0x1.52cc03d390dedp-29, -0x1.108a18f1099e3p-22, 0x1.a51d77be6c443p-16, // subdomain 4:7 (2.4375..2.5)
                -0x1.050e9aa28e5d0p-62, 0x1.f176f6d35dc08p-57, -0x1.cd7e248a61abep-50, 0x1.a06d704581fe0p-43, -0x1.6c7cdd774c2bbp-36, 0x1.34bb892083fdap-29, -0x1.f8a45eb138232p-23, 0x1.8c86e63446ae8p-16, // subdomain 4:8 (2.5..2.5625)
                -0x1.bd58123cbeecfp-63, 0x1.ad976a21577cdp-57, -0x1.93af5b670e18ep-50, 0x1.712edfa52efc2p-43, -0x1.47b62e39768eep-36, 0x1.19b8ebd7f7a68p-29, -0x1.d3c4c03f5b733p-23, 0x1.75bfcf9f135b3p-16, // subdomain 4:9 (2.5625..2.625)
                -0x1.7c9346aa22f3fp-63, 0x1.73a5245dbf9e0p-57, -0x1.61b85a480ee53p-50, 0x1.47d28117068d2p-43, -0x1.27150b1b4165ep-36, 0x1.016b7185bf275p-29, -0x1.b21946319a323p-23, 0x1.609fa4fe60506p-16, // subdomain 4:10 (2.625..2.6875)
                -0x1.45d37603de670p-63, 0x1.4214e44e88e42p-57, -0x1.36750dd01cbeap-50, 0x1.238d1c22c2e2ap-43, -0x1.0a1605f0061a8p-36, 0x1.d70ada84e295fp-30, -0x1.9350425695e36p-23, 0x1.4d01e86776513p-16, // subdomain 4:11 (2.6875..2.75)
                -0x1.17754d1c12731p-63, 0x1.179c5cd4fce2dp-57, -0x1.10ee946a6de27p-50, 0x1.03b0e090164c4p-43, -0x1.e08ea622ae8b0p-37, 0x1.af8583a25a632p-30, -0x1.77211f7317715p-23, 0x1.3ac5b85251b41p-16, // subdomain 4:12 (2.75..2.8125)
                -0x1.e03b0eabfa50ep-64, 0x1.e64dfa582ab19p-58, -0x1.e0a7b4b995248p-51, 0x1.cf51bec365366p-44, -0x1.b28c7f15e485cp-37, 0x1.8bcf3608cb69cp-30, -0x1.5d4b443dc1b3dp-23, 0x1.29cd6945013b3p-16, // subdomain 4:13 (2.8125..2.875)
                -0x1.9d57f34c36db3p-64, 0x1.a799801b714c4p-58, -0x1.a7e6f6174ff80p-51, 0x1.9dea80d954c27p-44, -0x1.897a15a1708c9p-37, 0x1.6b7eabcccb9edp-30, -0x1.45951c2d352dcp-23, 0x1.19fe2bf2722cdp-16, // subdomain 4:14 (2.875..2.9375)
                -0x1.646216a8a17aap-64, 0x1.7195400acf228p-58, -0x1.766d59e51b250p-51, 0x1.7250288eb965bp-44, -0x1.64c227c01495ap-37, 0x1.4e37210a037f1p-30, -0x1.2fcb40b591d69p-23, 0x1.0b3fbe2802c0ap-16, // subdomain 4:15 (2.9375..3.0)
                -0x1.33ca46d0977b7p-64, 0x1.42f8c4ab20bffp-58, -0x1.4b3a11951c9dap-51, 0x1.4bc55024900abp-44, -0x1.43e298d935810p-37, 0x1.33a6b6cbe6c16p-30, -0x1.1bbfbe6362529p-23, 0x1.faf84a419ad88p-17, // subdomain 4:16 (3.0..3.0625)
                -0x1.0a43542f2c48bp-64, 0x1.1aaef226acacfp-58, -0x1.25720a5a691b5p-51, 0x1.29a68aa051dc2p-44, -0x1.2669cdf479072p-37, 0x1.1b850f3b28a88p-30, -0x1.094971f190b2cp-23, 0x1.e13ee010761e2p-17, // subdomain 4:17 (3.0625..3.125)
                -0x1.cd6d438dfc22bp-65, 0x1.ef9b252bcd2a9p-59, -0x1.045a0e449186fp-51, 0x1.0b6697ae7fd2bp-44, -0x1.0bf46d4161b01p-37, 0x1.05921a99bc41fp-30, -0x1.f086f41770592p-24, 0x1.c92f0340edd32p-17, // subdomain 4:18 (3.125..3.1875)
                -0x1.90758daf7f6afp-65, 0x1.b31caeff4a739p-59, -0x1.cea3b726da680p-52, 0x1.e1165c192c548p-45, -0x1.e856e2fafb212p-38, 0x1.e32a1bd411c10p-31, -0x1.d11975d5438a8p-24, 0x1.b2a7c04c1eab0p-17, // subdomain 4:19 (3.1875..3.25)
                -0x1.5c17bcb435c6ep-65, 0x1.7e921abc8c70bp-59, -0x1.9ba00ff72eba2p-52, 0x1.b15487dbccc73p-45, -0x1.bd85089895a7dp-38, 0x1.beb6fa720ab3cp-31, -0x1.b40ee2c3a72b8p-24, 0x1.9d8b1552c08f0p-17, // subdomain 4:20 (3.25..3.3125)
                -0x1.2f0a68563d1d7p-65, 0x1.50dd951ab6dfbp-59, -0x1.6ebd1a50c7e23p-52, 0x1.86cf820276184p-45, -0x1.96ed295a38811p-38, 0x1.9d7128df95099p-31, -0x1.9931a69c008bcp-24, 0x1.89bda6625b2a9p-17, // subdomain 4:21 (3.3125..3.375)
                -0x1.08381e832466ep-65, 0x1.290bda4f3ad2ap-59, -0x1.472fc3e4ee30bp-52, 0x1.60e7125f88c0ep-45, -0x1.741940f20452ap-38, 0x1.7f08d4af89261p-31, -0x1.8051773ef06d8p-24, 0x1.77267a2189f48p-17, // subdomain 4:22 (3.375..3.4375)
                -0x1.cd6d92b8eced8p-66, 0x1.064d8719a1d0ep-59, -0x1.2449119163d56p-52, 0x1.3f0fa4087f4aap-45, -0x1.54a1528a9345ap-38, 0x1.6336d446ef2a5p-31, -0x1.6942c0f75cc50p-24, 0x1.65aebddf6068cp-17, // subdomain 4:23 (3.4375..3.5)
                -0x1.938124623a951p-66, 0x1.cfe3079b307f4p-60, -0x1.0571ffc79da2dp-52, 0x1.20cf6e00ade08p-45, -0x1.38299c7785b74p-38, 0x1.49bba0b4489a9p-31, -0x1.53de244e81915p-24, 0x1.55419023d7849p-17, // subdomain 4:24 (3.5..3.5625)
                -0x1.615bf9928bd13p-66, 0x1.9ac0b18c0ecb0p-60, -0x1.d4500ce2cc0c3p-53, 0x1.05bc062a3e63bp-45, -0x1.1e610af701a66p-38, 0x1.325e70f5cca49p-31, -0x1.4000033084813p-24, 0x1.45cbd0fc50ae2p-17, // subdomain 4:25 (3.5625..3.625)
                -0x1.35e2066d9c37ep-66, 0x1.6c3088c817ab2p-60, -0x1.a3f44b78ccbd7p-53, 0x1.daf09c34e4768p-46, -0x1.06ffe0d803284p-38, 0x1.1cec721f8c7b6p-31, -0x1.2d881b7009073p-24, 0x1.373bf75859dc3p-17, // subdomain 4:26 (3.625..3.6875)
                -0x1.1022096eda0d0p-66, 0x1.4353e3a92eb54p-60, -0x1.790ccb773eb44p-53, 0x1.af65582af91a6p-46, -0x1.e38d1c839d833p-39, 0x1.0938187b6954bp-31, -0x1.1c592cfac602cp-24, 0x1.2981eae0914c2p-17, // subdomain 4:27 (3.6875..3.75)
                -0x1.de9d8a4d2d24dp-67, 0x1.1f6b6b3fc3a12p-60, -0x1.52ef6d0c84d04p-53, 0x1.884709550c53cp-46, -0x1.bcf95e0e5720bp-39, 0x1.ee310c97cfa56p-32, -0x1.0c58aa4ae6feep-24, 0x1.1c8ee1c44b612p-17, // subdomain 4:28 (3.75..3.8125)
                -0x1.a572b29f78fccp-67, 0x1.ffa4e6faeabe5p-61, -0x1.31084318d227ap-53, 0x1.6517d67bdc9d2p-46, -0x1.99e055dc1fa3dp-39, 0x1.ccd20aa8c45c8p-32, -0x1.fadce3a518593p-25, 0x1.1055420adf019p-17, // subdomain 4:29 (3.8125..3.875)
                -0x1.7399133050d25p-67, 0x1.c7f61b3b91727p-61, -0x1.12d67ec4c621ep-53, 0x1.456902ae6acdfp-46, -0x1.79e8e654651f7p-39, 0x1.ae112160ac6a0p-32, -0x1.df091e928812fp-25, 0x1.04c88603973ffp-17, // subdomain 4:30 (3.875..3.9375)
                -0x1.481192c1602b2p-67, 0x1.96d56ca946190p-61, -0x1.efd39d6e225f4p-54, 0x1.28d8fb9d83c72p-46, -0x1.5cc3b1fa59361p-39, 0x1.91b2d9b6dc061p-32, -0x1.c50e07cfd7fe8p-25, 0x1.f3ba46f713ddfp-18, // subdomain 4:31 (3.9375..4.0)
                -0x1.10e06fed9be50p-58, 0x1.57c2e9b27b9e9p-53, -0x1.a98d91b6debd1p-47, 0x1.03212937568b7p-40, -0x1.35bd5e794b9c9p-34, 0x1.6b2a8af970fb2p-28, -0x1.a13ccbc444b2ep-22, 0x1.d5271052870b3p-16, // subdomain 5:0 (4.0..4.125)
                -0x1.ac946fe4b1a80p-59, 0x1.139b57e49ef4dp-53, -0x1.5c89dce66063ap-47, 0x1.b1dfa88497a98p-41, -0x1.093d73e2ced52p-34, 0x1.3e594c89ecc89p-28, -0x1.76b7161f8a3ecp-22, 0x1.b012d4ba44e8ap-16, // subdomain 5:1 (4.125..4.25)
                -0x1.522e872034f75p-59, 0x1.bbec1a3f27f8fp-54, -0x1.1ea7fc146b2bap-47, 0x1.6c9ecddf30bfap-41, -0x1.c7d52b7487481p-35, 0x1.17e7cc3fe622bp-28, -0x1.51628f83ac364p-22, 0x1.8ebb3c41f48fbp-16, // subdomain 5:2 (4.25..4.375)
                -0x1.0c1764d110ee3p-59, 0x1.671256c8c36b7p-54, -0x1.d96caf1b8d78fp-48, 0x1.338d40427c090p-41, -0x1.88fcb83f47d67p-35, 0x1.eda2d7a96f7c1p-29, -0x1.308411f33b52bp-22, 0x1.70acb1646e894p-16, // subdomain 5:3 (4.375..4.5)
                -0x1.aaf7644dee5c6p-60, 0x1.23a9356bf2fa7p-54, -0x1.8876422dda528p-48, 0x1.045667f206879p-41, -0x1.53e2c824b06eap-35, 0x1.b47edf2098eb1p-29, -0x1.137e12d55d9b2p-22, 0x1.5583bb6d69c7ep-16, // subdomain 5:4 (4.5..4.625)
                -0x1.5579a8ad5b63bp-60, 0x1.dbbcbfa8c97d8p-55, -0x1.4691416dae0fap-48, 0x1.ba40c1959228ap-42, -0x1.26dc764471ae6p-35, 0x1.83006c33dfc2dp-29, -0x1.f396f130f1cb0p-23, 0x1.3cea7815ecde2p-16, // subdomain 5:5 (4.625..4.75)
                -0x1.124050ff7d61bp-60, 0x1.8584f1938f0b1p-55, -0x1.10b8de7f01df8p-48, 0x1.78e224e0dd2c4p-42, -0x1.008e5f94c3285p-35, 0x1.58007e00de968p-29, -0x1.c5f6d67630074p-23, 0x1.269684e327cf0p-16, // subdomain 5:6 (4.75..4.875)
                -0x1.ba511a43a8973p-61, 0x1.40235c7e93660p-55, -0x1.c91b08c13d8d2p-49, 0x1.42330f11bb871p-42, -0x1.bfb964316752fp-36, 0x1.328934dffd88bp-29, -0x1.9d5bd909bfe81p-23, 0x1.1247434205207p-16, // subdomain 5:7 (4.875..5.0)
                -0x1.66183587b1382p-61, 0x1.0814f3275d268p-55, -0x1.805ed1bf8a74ap-49, 0x1.144bfa7d18181p-42, -0x1.87bd4a7144b1dp-36, 0x1.11cce846e993ep-29, -0x1.79221fbec772ep-23, 0x1.ff88cd3653072p-17, // subdomain 5:8 (5.0..5.125)
                -0x1.23043ca6b29a0p-61, 0x1.b53b146bdfd49p-56, -0x1.44441a3ee9f9fp-49, 0x1.db46213492bdbp-43, -0x1.57a92a8b28d64p-36, 0x1.ea3df5c24b692p-30, -0x1.58bd64a54148cp-23, 0x1.ddb97de2768e7p-17, // subdomain 5:9 (5.125..5.25)
                -0x1.dac1e1c5b1bf7p-62, 0x1.6b32eb8d9ef8ap-56, -0x1.126ce1e577debp-49, 0x1.99f1d8ba44672p-43, -0x1.2e4056f30bddbp-36, 0x1.b7dc1f5a9382dp-30, -0x1.3bb533954d42ap-23, 0x1.beca6a77f4d25p-17, // subdomain 5:10 (5.25..5.375)
                -0x1.84a47705e9195p-62, 0x1.2eb500bbeea89p-56, -0x1.d1e949164630fp-50, 0x1.62931c0f8893cp-43, -0x1.0a7c3e71b41f2p-36, 0x1.8b7eaf241f58ep-30, -0x1.21a1cfe5058b8p-23, 0x1.a26fe18d96552p-17, // subdomain 5:11 (5.375..5.5)
                -0x1.3f40e6cdf3854p-62, 0x1.fa355265ef275p-57, -0x1.8cacfab2b13c9p-50, 0x1.3381135b18505p-43, -0x1.d704e030daa42p-37, 0x1.64557d61683aap-30, -0x1.0a29a1ea5e4fcp-23, 0x1.88673d00025f7p-17, // subdomain 5:12 (5.5..5.625)
                -0x1.072326c6eb657p-62, 0x1.a894e2549fdd4p-57, -0x1.52b23394a640dp-50, 0x1.0b5f6a5778e6ep-43, -0x1.a1392ffd6a764p-37, 0x1.41af99d19b54dp-30, -0x1.e9fe266d25e72p-24, 0x1.7075a8eeb340fp-17, // subdomain 5:13 (5.625..5.75)
                -0x1.b33008e2a02f7p-63, 0x1.653223f036fc5p-57, -0x1.21fe3ab51e8b3p-50, 0x1.d21f1b8ef8565p-44, -0x1.726415bc62619p-37, 0x1.22f62cdb24a14p-30, -0x1.c3bd8ad9fa6f7p-24, 0x1.5a6719b012c15p-17, // subdomain 5:14 (5.75..5.875)
                -0x1.6901d9f06da8fp-63, 0x1.2d6327612a220p-57, -0x1.f1ebb5b7c3979p-51, 0x1.9749dfdd83463p-44, -0x1.49851c346bdedp-37, 0x1.07a8413b83d64p-30, -0x1.a11c27027b055p-24, 0x1.460d69078ffa5p-17, // subdomain 5:15 (5.875..6.0)
                -0x1.2c64a398900ffp-63, 0x1.fe0b7a4e4131dp-58, -0x1.ac9475371b53ap-51, 0x1.64b7581f4fcb9p-44, -0x1.25c43783d55d5p-37, 0x1.deae9411bf392p-31, -0x1.81b3a086ca8b7p-24, 0x1.333f9446b3d63p-17, // subdomain 5:16 (6.0..6.125)
                -0x1.f568089c1dc29p-64, 0x1.b0c5c1a4ec8c2p-58, -0x1.71d388c6cc244p-51, 0x1.39220e735880ap-44, -0x1.066adff6f48b6p-37, 0x1.b3488384153ebp-31, -0x1.652a6f33506b0p-24, 0x1.21d9161dd12dcp-17, // subdomain 5:17 (6.125..6.25)
                -0x1.a3ae09b1f7cf8p-64, 0x1.7030f1dc2ebd7p-58, -0x1.3fe86d29643c5p-51, 0x1.137ad3b483405p-44, -0x1.d5bce008a8401p-38, 0x1.8c7a8c615de32p-31, -0x1.4b321381b42a5p-24, 0x1.11b957c7635e0p-17, // subdomain 5:18 (6.25..6.375)
                -0x1.6043ae2e44460p-64, 0x1.3a0fe79e79010p-58, -0x1.1562538c87e21p-51, 0x1.e5be0f324e57ap-45, -0x1.a5370cdc96070p-38, 0x1.69b734da3c524p-31, -0x1.33859305524f9p-24, 0x1.02c335f51a5a5p-17, // subdomain 5:19 (6.375..6.5)
                -0x1.287ced1487111p-64, 0x1.0c91484fbba7ep-58, -0x1.e222118dedab1p-52, 0x1.ad2177075a6f1p-45, -0x1.7a647513361b4p-38, 0x1.4a83f96fe96bdp-31, -0x1.1de82f2e1b516p-24, 0x1.e9b92cfea3cf2p-18, // subdomain 5:20 (6.5..6.625)
                -0x1.f46a20f5ad167p-65, 0x1.cc74119f2974bp-59, -0x1.a3f1c7ddf994dp-52, 0x1.7be0a033aa31cp-45, -0x1.5485f5ab5e4d6p-38, 0x1.2e7682b7cd729p-31, -0x1.0a244cdc46e8dp-24, 0x1.cfdc18abc3bd3p-18, // subdomain 5:21 (6.625..6.75)
                -0x1.a764fbeec8e29p-65, 0x1.8ba8ecc06c7a4p-59, -0x1.6e9230c9e5da4p-52, 0x1.50ee8445a41f9p-45, -0x1.32f7ec2355e4dp-38, 0x1.15324d26ddfabp-31, -0x1.f01509d5766d3p-25, 0x1.b7c50f3235057p-18, // subdomain 5:22 (6.75..6.875)
                -0x1.672156ec1133ep-65, 0x1.54c586fa22339p-59, -0x1.40a7d7880a18bp-52, 0x1.2b67d610a823bp-45, -0x1.152dea8e72a55p-38, 0x1.fccd5ea4d70cfp-32, -0x1.cee1ae631cd94p-25, 0x1.a14e202f19ff6p-18, // subdomain 5:23 (6.875..7.0)
                -0x1.315dae21d0700p-65, 0x1.2628a9ded05e9p-59, -0x1.1910d82db8f99p-52, 0x1.0a8c0b6790166p-45, -0x1.f55e4903b0becp-39, 0x1.d39a5d530e62ep-32, -0x1.b063f53af292dp-25, 0x1.8c551168a2d28p-18, // subdomain 5:24 (7.0..7.125)
                -0x1.0444dd39cacdfp-65, 0x1.fcf4dfac2fb87p-60, -0x1.edb560d84bf5dp-53, 0x1.db6f48b6bb02ap-46, -0x1.c626e75f53ee7p-39, 0x1.ae502b1b0c38ap-32, -0x1.94598bc1e07c6p-25, 0x1.78baf4e90843dp-18, // subdomain 5:25 (7.125..7.25)
                -0x1.bcb151907d429p-66, 0x1.b93d2fc272cddp-60, -0x1.b273f384444c1p-53, 0x1.a8bee7f5ed98cp-46, -0x1.9c01b2ae93861p-39, 0x1.8c82800bbbf68p-32, -0x1.7a875a03b5e78p-25, 0x1.6663cc47c2512p-18, // subdomain 5:26 (7.25..7.375)
                -0x1.7cc17dee89298p-66, 0x1.7f5314d11e3c2p-60, -0x1.7f075851e6b6ep-53, 0x1.7c197511ebf74p-46, -0x1.7652cd788e8c8p-39, 0x1.6dd22c608b007p-32, -0x1.62b8a146e4a7dp-25, 0x1.55363752ed125p-18, // subdomain 5:27 (7.375..7.5)
                -0x1.46ba6d026c310p-66, 0x1.4daf011974c58p-60, -0x1.524f38e40e2bdp-53, 0x1.54b3d5256fc5cp-46, -0x1.5492f5c158a52p-39, 0x1.51eb5a30743f0p-32, -0x1.4cbe38da0245bp-25, 0x1.451b2c8f99e95p-18, // subdomain 5:28 (7.5..7.625)
                -0x1.18f7972c6bce9p-66, 0x1.230abca0dfc99p-60, -0x1.2b57f2aa5c2e9p-53, 0x1.31dff3a998373p-46, -0x1.364c8dffa07cdp-39, 0x1.38840fc53ca92p-32, -0x1.386de4b25bcfdp-25, 0x1.35fdba418e3c5p-18, // subdomain 5:29 (7.625..7.75)
                -0x1.e43d15d4c0ccbp-67, 0x1.fcab9f4ad9d8ep-61, -0x1.09532d262c70ap-53, 0x1.130849a254606p-46, -0x1.1b191bccb85f7p-39, 0x1.215ae88265931p-32, -0x1.25a1c22456ad5p-25, 0x1.27cacedb8feb9p-18, // subdomain 5:30 (7.75..7.875)
                -0x1.a222426689c5ep-67, 0x1.bd58317f05824p-61, -0x1.d7237d44932d2p-54, 0x1.ef5842fcbfbe7p-47, -0x1.029f28c6e928dp-39, 0x1.0c35fbf2f047bp-32, -0x1.1437c791f2b35p-25, 0x1.1a7107e26666dp-18, // subdomain 5:31 (7.875..8.0)
                -0x1.50e37c46b75fcp-58, 0x1.6e380eaa6ffeap-53, -0x1.8b386706d8352p-47, 0x1.a87597e9e5a94p-41, -0x1.c4b720d629908p-35, 0x1.dfaf7e5de5748p-29, -0x1.f8dd638965aa9p-23, 0x1.07df2a8541e55p-16, // subdomain 6:0 (8.0..8.25)
                -0x1.fc3690c76d3eap-59, 0x1.1bb2a4e231abcp-53, -0x1.3a8e10f029dd9p-47, 0x1.5b2ade598985cp-41, -0x1.7ca69e3f77fb0p-35, 0x1.9ec46e723e766p-29, -0x1.c11425bd3bdb1p-23, 0x1.e3173af5b3ec6p-17, // subdomain 6:1 (8.25..8.5)
                -0x1.8220461aca39cp-59, 0x1.ba7944e46899fp-54, -0x1.f7b91ffb9e9b1p-48, 0x1.1d790ddcdc36ep-41, -0x1.418e122a44dbcp-35, 0x1.680ed6f1c92e9p-29, -0x1.90bf8c6587428p-23, 0x1.bb4efa1509daap-17, // subdomain 6:2 (8.5..8.75)
                -0x1.27677f3c7ddd0p-59, 0x1.5b3ed522599c5p-54, -0x1.95a33bbe72c92p-48, 0x1.d7e069d0263c4p-42, -0x1.10d7169ae7820p-35, 0x1.39bdde89f6a30p-29, -0x1.66b709b950186p-23, 0x1.97c0fdbde918cp-17, // subdomain 6:3 (8.75..9.0)
                -0x1.c6fcc9560858ap-60, 0x1.1228eaac415c2p-54, -0x1.486ee158bb9e3p-48, 0x1.87e3809dc020dp-42, -0x1.d0f752bdf2c7fp-36, 0x1.125c6164eb019p-29, -0x1.420633ad4774dp-23, 0x1.77e1f5fcd0730p-17, // subdomain 6:4 (9.0..9.25)
                -0x1.609a789caea04p-60, 0x1.b368f1ca96560p-55, -0x1.0b4e2224137b0p-48, 0x1.46f65e86e3151p-42, -0x1.8dc85ab5d9396p-36, 0x1.e17a5570af075p-30, -0x1.21e297cfbbd3bp-23, 0x1.5b3b733d35b82p-17, // subdomain 6:5 (9.25..9.5)
                -0x1.12e83108b8cd5p-60, 0x1.5ba75e512a7b7p-55, -0x1.b5443cf05763dp-49, 0x1.11fe5f832f12cp-42, -0x1.559bf98334661p-36, 0x1.a7d793a2912dbp-30, -0x1.05a3baaacd048p-23, 0x1.41685ae40b3ccp-17, // subdomain 6:6 (9.5..9.75)
                -0x1.af2267bcdcc91p-61, 0x1.170b6e8a98289p-55, -0x1.67574c0495f09p-49, 0x1.cd241c09a394ep-43, -0x1.267082673abecp-36, 0x1.7641b8c9945f7p-30, -0x1.d9798529b645ap-24, 0x1.2a120712e6acbp-17, // subdomain 6:7 (9.75..10.0)
                -0x1.53f017a94b58ep-61, 0x1.c23430871037cp-56, -0x1.28a39bd219e0dp-49, 0x1.859d93779c591p-43, -0x1.fd55b2a8ca395p-37, 0x1.4b71ebd6293e2p-30, -0x1.ad6edaf2b5e38p-24, 0x1.14edee3cf140dp-17, // subdomain 6:8 (10.0..10.25)
                -0x1.0d72bd1cc48d8p-61, 0x1.6cec126199346p-56, -0x1.ebe06911e476bp-50, 0x1.4a72662a17116p-43, -0x1.ba008ae257f4bp-37, 0x1.265b07e693a65p-30, -0x1.866010f0fc012p-24, 0x1.01bbb71bdef16p-17, // subdomain 6:9 (10.25..10.5)
                -0x1.ad50528ff26acp-62, 0x1.2928bdb05be3fp-56, -0x1.997f26689fedfp-50, 0x1.194aecae09ec2p-43, -0x1.80cad18273b64p-37, 0x1.061e70394c9e8p-30, -0x1.63a463201748cp-24, 0x1.e08747286b530p-18, // subdomain 6:10 (10.5..10.75)
                -0x1.57acee0e9768dp-62, 0x1.e618f53495e22p-57, -0x1.5644926a93341p-50, 0x1.e0958b216a081p-44, -0x1.50022f2d112d0p-37, 0x1.d40678a693d2ap-31, -0x1.44ac70ea40bdbp-24, 0x1.c0aa83610bf6ap-18, // subdomain 6:11 (10.75..11.0)
                -0x1.146714d2917a8p-62, 0x1.8f44d5196ca63p-57, -0x1.1f2a78d994778p-50, 0x1.9beb373ec28d6p-44, -0x1.26443eede2dabp-37, 0x1.a2de6a3221171p-31, -0x1.28fdfb0f36149p-24, 0x1.a38caac9f1144p-18, // subdomain 6:12 (11.0..11.25)
                -0x1.be9648694d802p-63, 0x1.49490091b207dp-57, -0x1.e3a193bfc3f28p-51, 0x1.6234aa7ac5a3ep-44, -0x1.026e95917b64fp-37, 0x1.77c278288b3fcp-31, -0x1.10306a56d0b67p-24, 0x1.88e3e114140fbp-18, // subdomain 6:13 (11.25..11.5)
                -0x1.6a54a861954cap-63, 0x1.10a169eba9f52p-57, -0x1.98af2111a8842p-51, 0x1.318626f350007p-44, -0x1.c7246e5fafc19p-38, 0x1.51d99caff1621p-31, -0x1.f3d3edd4888abp-25, 0x1.706f613f2be8ep-18, // subdomain 6:14 (11.5..11.75)
                -0x1.27300b28f2c86p-63, 0x1.c5251f5a91f97p-58, -0x1.5a8448cf5e8b8p-51, 0x1.0851a643db73fp-44, -0x1.91d363544c6c4p-38, 0x1.306ccd44bf77ap-31, -0x1.cbbaa1a03f7cap-25, 0x1.59f6372a02e34p-18, // subdomain 6:15 (11.75..12.0)
            };

            alignas(64) const double coeffs1[208 * 2] = {
                -0x1.d07126526127fp-10, 0x1.be2e5581c20e9p-1, // subdomain 0:0 (0.125..0.12890625)
                -0x1.cd6df235069cfp-10, 0x1.bc5f667f3c6dbp-1, // subdomain 0:1 (0.12890625..0.1328125)
                -0x1.ca712ad4af8c7p-10, 0x1.ba93777919864p-1, // subdomain 0:2 (0.1328125..0.13671875)
                -0x1.c77abff758f19p-10, 0x1.b8ca820abd3d7p-1, // subdomain 0:3 (0.13671875..0.140625)
                -0x1.c48aa1918ac3cp-10, 0x1.b7047fdfae29ap-1, // subdomain 0:4 (0.140625..0.14453125)
                -0x1.c1a0bfc5c428ep-10, 0x1.b5416ab367513p-1, // subdomain 0:5 (0.14453125..0.1484375)
                -0x1.bebd0ae3e9a62p-10, 0x1.b3813c512a7cdp-1, // subdomain 0:6 (0.1484375..0.15234375)
                -0x1.bbdf7368b548bp-10, 0x1.b1c3ee93d319dp-1, // subdomain 0:7 (0.15234375..0.15625)
                -0x1.b907e9fd28be3p-10, 0x1.b0097b65a9ac0p-1, // subdomain 0:8 (0.15625..0.16015625)
                -0x1.b6365f7601567p-10, 0x1.ae51dcc037cc2p-1, // subdomain 0:9 (0.16015625..0.1640625)
                -0x1.b36ac4d32de6cp-10, 0x1.ac9d0cac1cb16p-1, // subdomain 0:10 (0.1640625..0.16796875)
                -0x1.b0a50b3f46888p-10, 0x1.aaeb0540e2442p-1, // subdomain 0:11 (0.16796875..0.171875)
                -0x1.ade5240f062b9p-10, 0x1.a93bc0a4d2b85p-1, // subdomain 0:12 (0.171875..0.17578125)
                -0x1.ab2b00c0c5f64p-10, 0x1.a78f390cceacbp-1, // subdomain 0:13 (0.17578125..0.1796875)
                -0x1.a87692fbfa6bbp-10, 0x1.a5e568bc23ce6p-1, // subdomain 0:14 (0.1796875..0.18359375)
                -0x1.a5c7cc90b2525p-10, 0x1.a43e4a0463fe2p-1, // subdomain 0:15 (0.18359375..0.1875)
                -0x1.a31e9f7717547p-10, 0x1.a299d7453cf60p-1, // subdomain 0:16 (0.1875..0.19140625)
                -0x1.a07afdcef0535p-10, 0x1.a0f80aec506d7p-1, // subdomain 0:17 (0.19140625..0.1953125)
                -0x1.9ddcd9df25684p-10, 0x1.9f58df750cb9fp-1, // subdomain 0:18 (0.1953125..0.19921875)
                -0x1.9b442615458bep-10, 0x1.9dbc4f6885eb8p-1, // subdomain 0:19 (0.19921875..0.203125)
                -0x1.98b0d5050dde7p-10, 0x1.9c22555d4f621p-1, // subdomain 0:20 (0.203125..0.20703125)
                -0x1.9622d967f28b2p-10, 0x1.9a8aebf755dafp-1, // subdomain 0:21 (0.20703125..0.2109375)
                -0x1.939a261ca9405p-10, 0x1.98f60de7b9f50p-1, // subdomain 0:22 (0.2109375..0.21484375)
                -0x1.9116ae26b536fp-10, 0x1.9763b5ecab28dp-1, // subdomain 0:23 (0.21484375..0.21875)
                -0x1.8e9864adf4c35p-10, 0x1.95d3ded14334dp-1, // subdomain 0:24 (0.21875..0.22265625)
                -0x1.8c1f3cfe3069fp-10, 0x1.9446836d61fadp-1, // subdomain 0:25 (0.22265625..0.2265625)
                -0x1.89ab2a86ab72dp-10, 0x1.92bb9ea589cd4p-1, // subdomain 0:26 (0.2265625..0.23046875)
                -0x1.873c20d9b5f54p-10, 0x1.91332b6abc2bfp-1, // subdomain 0:27 (0.23046875..0.234375)
                -0x1.84d213ac40577p-10, 0x1.8fad24ba56ed5p-1, // subdomain 0:28 (0.234375..0.23828125)
                -0x1.826cf6d5703cbp-10, 0x1.8e29859df1d43p-1, // subdomain 0:29 (0.23828125..0.2421875)
                -0x1.800cbe4e36dbep-10, 0x1.8ca8492b3c8efp-1, // subdomain 0:30 (0.2421875..0.24609375)
                -0x1.7db15e30e8ba8p-10, 0x1.8b296a83dd209p-1, // subdomain 0:31 (0.24609375..0.25)
                -0x1.7a314a1639a80p-9, 0x1.88ef81dcc0893p-1, // subdomain 1:0 (0.25..0.2578125)
                -0x1.75970d49d68d0p-9, 0x1.85ffbca53afd5p-1, // subdomain 1:1 (0.2578125..0.265625)
                -0x1.710f6268654c5p-9, 0x1.831919465cfecp-1, // subdomain 1:2 (0.265625..0.2734375)
                -0x1.6c99f043f98d6p-9, 0x1.803b72f5d0e73p-1, // subdomain 1:3 (0.2734375..0.28125)
                -0x1.68365f96ab609p-9, 0x1.7d66a599b2f06p-1, // subdomain 1:4 (0.28125..0.2890625)
                -0x1.63e45af7036e5p-9, 0x1.7a9a8dc4cccbap-1, // subdomain 1:5 (0.2890625..0.296875)
                -0x1.5fa38eccb3d13p-9, 0x1.77d708b2e814ap-1, // subdomain 1:6 (0.296875..0.3046875)
                -0x1.5b73a9459b770p-9, 0x1.751bf44537078p-1, // subdomain 1:7 (0.3046875..0.3125)
                -0x1.57545a4b11f78p-9, 0x1.72692efed2e6ep-1, // subdomain 1:8 (0.3125..0.3203125)
                -0x1.5345537779ddep-9, 0x1.6fbe98014f844p-1, // subdomain 1:9 (0.3203125..0.328125)
                -0x1.4f46480c17711p-9, 0x1.6d1c0f09635e3p-1, // subdomain 1:10 (0.328125..0.3359375)
                -0x1.4b56ece72a171p-9, 0x1.6a81746ba3cedp-1, // subdomain 1:11 (0.3359375..0.34375)
                -0x1.4776f87a467c4p-9, 0x1.67eea91154c59p-1, // subdomain 1:12 (0.34375..0.3515625)
                -0x1.43a622c0efb67p-9, 0x1.65638e754b8e6p-1, // subdomain 1:13 (0.3515625..0.359375)
                -0x1.3fe425376da87p-9, 0x1.62e006a0e4283p-1, // subdomain 1:14 (0.359375..0.3671875)
                -0x1.3c30bad1def89p-9, 0x1.6063f42908b34p-1, // subdomain 1:15 (0.3671875..0.375)
                -0x1.388b9ff384f7ap-9, 0x1.5def3a2b4a816p-1, // subdomain 1:16 (0.375..0.3828125)
                -0x1.34f4926647e3bp-9, 0x1.5b81bc4b0c554p-1, // subdomain 1:17 (0.3828125..0.390625)
                -0x1.316b515271fd9p-9, 0x1.591b5eaebd630p-1, // subdomain 1:18 (0.390625..0.3984375)
                -0x1.2def9d369ff1bp-9, 0x1.56bc05fd24a40p-1, // subdomain 1:19 (0.3984375..0.40625)
                -0x1.2a8137dfe5230p-9, 0x1.5463975abc16cp-1, // subdomain 1:20 (0.40625..0.4140625)
                -0x1.271fe462226e8p-9, 0x1.5211f8671b829p-1, // subdomain 1:21 (0.4140625..0.421875)
                -0x1.23cb67108e0a4p-9, 0x1.4fc70f3a725d4p-1, // subdomain 1:22 (0.421875..0.4296875)
                -0x1.208385766b2c7p-9, 0x1.4d82c26310702p-1, // subdomain 1:23 (0.4296875..0.4375)
                -0x1.1d48064ff02ffp-9, 0x1.4b44f8e2fcdffp-1, // subdomain 1:24 (0.4375..0.4453125)
                -0x1.1a18b18359f72p-9, 0x1.490d9a2d9b3a7p-1, // subdomain 1:25 (0.4453125..0.453125)
                -0x1.16f5501a2b54cp-9, 0x1.46dc8e255e312p-1, // subdomain 1:26 (0.453125..0.4609375)
                -0x1.13ddac3a974d5p-9, 0x1.44b1bd1987aa4p-1, // subdomain 1:27 (0.4609375..0.46875)
                -0x1.10d19121150bdp-9, 0x1.428d0fc3f5d36p-1, // subdomain 1:28 (0.46875..0.4765625)
                -0x1.0dd0cb1a1c6bap-9, 0x1.406e6f46fce2ep-1, // subdomain 1:29 (0.4765625..0.484375)
                -0x1.0adb277c0a038p-9, 0x1.3e55c52b4d398p-1, // subdomain 1:30 (0.484375..0.4921875)
                -0x1.07f074a129a47p-9, 0x1.3c42fb5de5946p-1, // subdomain 1:31 (0.4921875..0.5)
                -0x1.03a4819c57725p-8, 0x1.3931a1e46a521p-1, // subdomain 2:0 (0.5..0.515625)
                -0x1.fc1ce9c776aabp-9, 0x1.352e49b054268p-1, // subdomain 2:1 (0.515625..0.53125)
                -0x1.f141bcaed92d2p-9, 0x1.3140f84b0a0e0p-1, // subdomain 2:2 (0.53125..0.546875)
                -0x1.e6b4bb9bbedbcp-9, 0x1.2d690ea0dbf16p-1, // subdomain 2:3 (0.546875..0.5625)
                -0x1.dc7341b19e96ap-9, 0x1.29a5f3031528ap-1, // subdomain 2:4 (0.5625..0.578125)
                -0x1.d27ac458d8a75p-9, 0x1.25f710f243800p-1, // subdomain 2:5 (0.578125..0.59375)
                -0x1.c8c8d21b9b795p-9, 0x1.225bd8ead2dc4p-1, // subdomain 2:6 (0.59375..0.609375)
                -0x1.bf5b119092f5ep-9, 0x1.1ed3c033e12c5p-1, // subdomain 2:7 (0.609375..0.625)
                -0x1.b62f4052ade75p-9, 0x1.1b5e40b02fc8bp-1, // subdomain 2:8 (0.625..0.640625)
                -0x1.ad4332054ea6ap-9, 0x1.17fad8b118c10p-1, // subdomain 2:9 (0.640625..0.65625)
                -0x1.a494cf6445a3cp-9, 0x1.14a90acb6fe4bp-1, // subdomain 2:10 (0.65625..0.671875)
                -0x1.9c22155efd1ebp-9, 0x1.11685dae388c0p-1, // subdomain 2:11 (0.671875..0.6875)
                -0x1.93e9143e45a7cp-9, 0x1.0e385bfb1a51dp-1, // subdomain 2:12 (0.6875..0.703125)
                -0x1.8be7eed43ac99p-9, 0x1.0b18942080064p-1, // subdomain 2:13 (0.703125..0.71875)
                -0x1.841cd9b5be948p-9, 0x1.080898354d2bap-1, // subdomain 2:14 (0.71875..0.734375)
                -0x1.7c861a7d11b04p-9, 0x1.0507fdd6174c8p-1, // subdomain 2:15 (0.734375..0.75)
                -0x1.75220715141f5p-9, 0x1.02165e03d165ap-1, // subdomain 2:16 (0.75..0.765625)
                -0x1.6def050cc1093p-9, 0x1.fe66aa07b103fp-2, // subdomain 2:17 (0.765625..0.78125)
                -0x1.66eb88f27daf5p-9, 0x1.f8bd0482a2f5ep-2, // subdomain 2:18 (0.78125..0.796875)
                -0x1.601615b6d92c0p-9, 0x1.f32f105f9133cp-2, // subdomain 2:19 (0.796875..0.8125)
                -0x1.596d3c165fca0p-9, 0x1.edbc185e092f1p-2, // subdomain 2:20 (0.8125..0.828125)
                -0x1.52ef9a0a29a2cp-9, 0x1.e8636cdffe567p-2, // subdomain 2:21 (0.828125..0.84375)
                -0x1.4c9bda3ed0d21p-9, 0x1.e32463b5ec0edp-2, // subdomain 2:22 (0.84375..0.859375)
                -0x1.4670b3917fe85p-9, 0x1.ddfe57ed0f89dp-2, // subdomain 2:23 (0.859375..0.875)
                -0x1.406ce892cd605p-9, 0x1.d8f0a99fa1b02p-2, // subdomain 2:24 (0.875..0.890625)
                -0x1.3a8f470f1cc60p-9, 0x1.d3fabdc6fa7b2p-2, // subdomain 2:25 (0.890625..0.90625)
                -0x1.34d6a79c41e38p-9, 0x1.cf1bfe0f88300p-2, // subdomain 2:26 (0.90625..0.921875)
                -0x1.2f41ed2c25c67p-9, 0x1.ca53d8ae85f5bp-2, // subdomain 2:27 (0.921875..0.9375)
                -0x1.29d004a430b74p-9, 0x1.c5a1c0395e392p-2, // subdomain 2:28 (0.9375..0.953125)
                -0x1.247fe4793f578p-9, 0x1.c1052b7ea6406p-2, // subdomain 2:29 (0.953125..0.96875)
                -0x1.1f508c4fec0b3p-9, 0x1.bc7d9560a12f2p-2, // subdomain 2:30 (0.96875..0.984375)
                -0x1.1a4104a0f8998p-9, 0x1.b80a7cb139941p-2, // subdomain 2:31 (0.984375..1.0)
                -0x1.12e35755d082ep-8, 0x1.b18331970ef07p-2, // subdomain 3:0 (1.0..1.03125)
                -0x1.0978011ad5cc5p-8, 0x1.a91210b28a0ecp-2, // subdomain 3:1 (1.03125..1.0625)
                -0x1.007c1e0967a43p-8, 0x1.a0ea885c67852p-2, // subdomain 3:2 (1.0625..1.09375)
                -0x1.efd2c6911df50p-9, 0x1.9909366fd1126p-2, // subdomain 3:3 (1.09375..1.125)
                -0x1.df73df0356426p-9, 0x1.916ae97347c55p-2, // subdomain 3:4 (1.125..1.15625)
                -0x1.cfd08636a4ffbp-9, 0x1.8a0c9d6c1656dp-2, // subdomain 3:5 (1.15625..1.1875)
                -0x1.c0de722afc48fp-9, 0x1.82eb78ed5d634p-2, // subdomain 3:6 (1.1875..1.21875)
                -0x1.b29400f19fee3p-9, 0x1.7c04ca5de8de1p-2, // subdomain 3:7 (1.21875..1.25)
                -0x1.a4e82c983448ap-9, 0x1.755605706a28bp-2, // subdomain 3:8 (1.25..1.28125)
                -0x1.97d28009056eep-9, 0x1.6edcc0ca141a2p-2, // subdomain 3:9 (1.28125..1.3125)
                -0x1.8b4b0ccb67f66p-9, 0x1.6896b3d3ef031p-2, // subdomain 3:10 (1.3125..1.34375)
                -0x1.7f4a619108a5dp-9, 0x1.6281b4b38b3b0p-2, // subdomain 3:11 (1.34375..1.375)
                -0x1.73c9817ec369fp-9, 0x1.5c9bb66801da0p-2, // subdomain 3:12 (1.375..1.40625)
                -0x1.68c1dc2130a6cp-9, 0x1.56e2c70875a32p-2, // subdomain 3:13 (1.40625..1.4375)
                -0x1.5e2d45fe87370p-9, 0x1.51550e21827e9p-2, // subdomain 3:14 (1.4375..1.46875)
                -0x1.5405f1b8bd226p-9, 0x1.4bf0cb2f40c57p-2, // subdomain 3:15 (1.46875..1.5)
                -0x1.4a4669b3fcedap-9, 0x1.46b45431b397ep-2, // subdomain 3:16 (1.5..1.53125)
                -0x1.40e98a3694a13p-9, 0x1.419e1459a6f68p-2, // subdomain 3:17 (1.53125..1.5625)
                -0x1.37ea7bf87920cp-9, 0x1.3cac8acc2bd0cp-2, // subdomain 3:18 (1.5625..1.59375)
                -0x1.2f44af1956c7ep-9, 0x1.37de497b060c4p-2, // subdomain 3:19 (1.59375..1.625)
                -0x1.26f3d674f1c3ap-9, 0x1.3331f410830a3p-2, // subdomain 3:20 (1.625..1.65625)
                -0x1.1ef3e34e4f610p-9, 0x1.2ea63eed4eb44p-2, // subdomain 3:21 (1.65625..1.6875)
                -0x1.1741014ac75dfp-9, 0x1.2a39ee36f9e2cp-2, // subdomain 3:22 (1.6875..1.71875)
                -0x1.0fd792b6b2e23p-9, 0x1.25ebd4f5ff413p-2, // subdomain 3:23 (1.71875..1.75)
                -0x1.08b42d0df7965p-9, 0x1.21bad4422be96p-2, // subdomain 3:24 (1.75..1.78125)
                -0x1.01d395c329b39p-9, 0x1.1da5da7c66f55p-2, // subdomain 3:25 (1.78125..1.8125)
                -0x1.f6657e80e2216p-10, 0x1.19abe294e7792p-2, // subdomain 3:26 (1.8125..1.84375)
                -0x1.e99d8c3b8864ep-10, 0x1.15cbf35cfad36p-2, // subdomain 3:27 (1.84375..1.875)
                -0x1.dd49dd10d46b0p-10, 0x1.12051ee38e372p-2, // subdomain 3:28 (1.875..1.90625)
                -0x1.d16543b013d70p-10, 0x1.0e5681dbbde45p-2, // subdomain 3:29 (1.90625..1.9375)
                -0x1.c5ead89bb4f4ep-10, 0x1.0abf430cbac83p-2, // subdomain 3:30 (1.9375..1.96875)
                -0x1.bad5f5f5a1e56p-10, 0x1.073e92ca64614p-2, // subdomain 3:31 (1.96875..2.0)
                -0x1.aaeb6de6cd87ap-9, 0x1.0226258f7ee2dp-2, // subdomain 4:0 (2.0..2.0625)
                -0x1.96eed1c6ff575p-9, 0x1.f745ca538915bp-3, // subdomain 4:1 (2.0625..2.125)
                -0x1.84400062608bdp-9, 0x1.ead9e428d6984p-3, // subdomain 4:2 (2.125..2.1875)
                -0x1.72c3c2647ea03p-9, 0x1.defe98ffc98e1p-3, // subdomain 4:3 (2.1875..2.25)
                -0x1.62617ebe608b9p-9, 0x1.d3aab84699bd9p-3, // subdomain 4:4 (2.25..2.3125)
                -0x1.5302f136393eap-9, 0x1.c8d5ccb2d0723p-3, // subdomain 4:5 (2.3125..2.375)
                -0x1.4493e9ea28bd0p-9, 0x1.be780aa21fdd2p-3, // subdomain 4:6 (2.375..2.4375)
                -0x1.37021494b8ad4p-9, 0x1.b48a405f617e6p-3, // subdomain 4:7 (2.4375..2.5)
                -0x1.2a3cc68db2093p-9, 0x1.ab05c811de297p-3, // subdomain 4:8 (2.5..2.5625)
                -0x1.1e34d2a71729cp-9, 0x1.a1e47b2494758p-3, // subdomain 4:9 (2.5625..2.625)
                -0x1.12dc622592727p-9, 0x1.9920a6f9b28a1p-3, // subdomain 4:10 (2.625..2.6875)
                -0x1.0826d22e707fap-9, 0x1.90b502c40fb7bp-3, // subdomain 4:11 (2.6875..2.75)
                -0x1.fc112a380e0ffp-10, 0x1.889ca66543fd9p-3, // subdomain 4:12 (2.75..2.8125)
                -0x1.e8ee2e79b21c7p-10, 0x1.80d3023324c7ap-3, // subdomain 4:13 (2.8125..2.875)
                -0x1.d6d14d24c7eacp-10, 0x1.7953d78b07863p-3, // subdomain 4:14 (2.875..2.9375)
                -0x1.c5a8ba656f761p-10, 0x1.721b321c4911bp-3, // subdomain 4:15 (2.9375..3.0)
                -0x1.b5641bf859ce5p-10, 0x1.6b2561d64dbbep-3, // subdomain 4:16 (3.0..3.0625)
                -0x1.a5f466aa0682bp-10, 0x1.646ef568886afp-3, // subdomain 4:17 (3.0625..3.125)
                -0x1.974bbf736e360p-10, 0x1.5df4b54523d3dp-3, // subdomain 4:18 (3.125..3.1875)
                -0x1.895d5fc93f070p-10, 0x1.57b39f18b2a25p-3, // subdomain 4:19 (3.1875..3.25)
                -0x1.7c1d7cc039a54p-10, 0x1.51a8e1aadb528p-3, // subdomain 4:20 (3.25..3.3125)
                -0x1.6f8130b4d8460p-10, 0x1.4bd1d91d527d3p-3, // subdomain 4:21 (3.3125..3.375)
                -0x1.637e672e76076p-10, 0x1.462c0b7fa9219p-3, // subdomain 4:22 (3.375..3.4375)
                -0x1.580bcabee64fep-10, 0x1.40b525af81bfdp-3, // subdomain 4:23 (3.4375..3.5)
                -0x1.4d20b4a7006aep-10, 0x1.3b6af87dabd6ap-3, // subdomain 4:24 (3.5..3.5625)
                -0x1.42b51e0f3cea8p-10, 0x1.364b761175c59p-3, // subdomain 4:25 (3.5625..3.625)
                -0x1.38c192a942fdap-10, 0x1.3154af843cd99p-3, // subdomain 4:26 (3.625..3.6875)
                -0x1.2f3f249449a68p-10, 0x1.2c84d2afe58d9p-3, // subdomain 4:27 (3.6875..3.75)
                -0x1.2627616294eb1p-10, 0x1.27da282b757c4p-3, // subdomain 4:28 (3.75..3.8125)
                -0x1.1d7448222f70ep-10, 0x1.2353117187869p-3, // subdomain 4:29 (3.8125..3.875)
                -0x1.1520404e60166p-10, 0x1.1eee072cc349ep-3, // subdomain 4:30 (3.875..3.9375)
                -0x1.0d26119250fbfp-10, 0x1.1aa997a6e4f8ep-3, // subdomain 4:31 (3.9375..4.0)
                -0x1.01ccafd71bfc1p-9, 0x1.147d1b190ca46p-3, // subdomain 5:0 (4.0..4.125)
                -0x1.e774e0c5f564dp-10, 0x1.0ca7c44d709a0p-3, // subdomain 5:1 (4.125..4.25)
                -0x1.cd831a6b740b8p-10, 0x1.053e861ffc32bp-3, // subdomain 5:2 (4.25..4.375)
                -0x1.b58bf67a385f3p-10, 0x1.fc721085dd723p-4, // subdomain 5:3 (4.375..4.5)
                -0x1.9f5e138ccd942p-10, 0x1.ef1f8a0371197p-4, // subdomain 5:4 (4.5..4.625)
                -0x1.8acdd807fd07ep-10, 0x1.e277e0aabb532p-4, // subdomain 5:5 (4.625..4.75)
                -0x1.77b4a9eed63d8p-10, 0x1.d66ec4c18c3efp-4, // subdomain 5:6 (4.75..4.875)
                -0x1.65f04531f3666p-10, 0x1.caf9099dc8c46p-4, // subdomain 5:7 (4.875..5.0)
                -0x1.55622b5951f3bp-10, 0x1.c00c8546be50ap-4, // subdomain 5:8 (5.0..5.125)
                -0x1.45ef285897fc9p-10, 0x1.b59ff441f8d00p-4, // subdomain 5:9 (5.125..5.25)
                -0x1.377ee91ea83d4p-10, 0x1.abaae0f07717ep-4, // subdomain 5:10 (5.25..5.375)
                -0x1.29fba10cf878fp-10, 0x1.a2258dfa0e771p-4, // subdomain 5:11 (5.375..5.5)
                -0x1.1d51bc00fe192p-10, 0x1.9908e359616dfp-4, // subdomain 5:12 (5.5..5.625)
                -0x1.116f9b0042596p-10, 0x1.904e5dabd3f12p-4, // subdomain 5:13 (5.625..5.75)
                -0x1.064559eb0983ep-10, 0x1.87efff770a69ap-4, // subdomain 5:14 (5.75..5.875)
                -0x1.f78939bb5783cp-11, 0x1.7fe8442147739p-4, // subdomain 5:15 (5.875..6.0)
                -0x1.e3c0c8417a15cp-11, 0x1.78321463d4233p-4, // subdomain 5:16 (6.0..6.125)
                -0x1.d119cb6c875d5p-11, 0x1.70c8bc06df3b9p-4, // subdomain 5:17 (6.125..6.25)
                -0x1.bf7ed770bd46fp-11, 0x1.69a7e0bd3074fp-4, // subdomain 5:18 (6.25..6.375)
                -0x1.aedc6dbab6490p-11, 0x1.62cb79fbe7819p-4, // subdomain 5:19 (6.375..6.5)
                -0x1.9f20c91f0bab0p-11, 0x1.5c2fc9af72218p-4, // subdomain 5:20 (6.5..6.625)
                -0x1.903bb039e8f5dp-11, 0x1.55d155b317cc7p-4, // subdomain 5:21 (6.625..6.75)
                -0x1.821e4d2bc1ef9p-11, 0x1.4face1f40bfa1p-4, // subdomain 5:22 (6.75..6.875)
                -0x1.74bb09ff2b3a8p-11, 0x1.49bf6b2c03e4cp-4, // subdomain 5:23 (6.875..7.0)
                -0x1.6805711dc16efp-11, 0x1.44062221e74f4p-4, // subdomain 5:24 (7.0..7.125)
                -0x1.5bf2114fe8b15p-11, 0x1.3e7e67616f305p-4, // subdomain 5:25 (7.125..7.25)
                -0x1.507664d5ade2dp-11, 0x1.3925c75c6ffa9p-4, // subdomain 5:26 (7.25..7.375)
                -0x1.4588bb36c1c40p-11, 0x1.33f9f6ea33d88p-4, // subdomain 5:27 (7.375..7.5)
                -0x1.3b202574db8d1p-11, 0x1.2ef8d01ab5f30p-4, // subdomain 5:28 (7.5..7.625)
                -0x1.313464572cd5ap-11, 0x1.2a204f54cb4e3p-4, // subdomain 5:29 (7.625..7.75)
                -0x1.27bdd88f4494dp-11, 0x1.256e90b757390p-4, // subdomain 5:30 (7.75..7.875)
                -0x1.1eb574800d4adp-11, 0x1.20e1cdb696c9dp-4, // subdomain 5:31 (7.875..8.0)
                -0x1.11e93f59f33bap-10, 0x1.1a50610170102p-4, // subdomain 6:0 (8.0..8.25)
                -0x1.0221654dd3f7cp-10, 0x1.12012483654f3p-4, // subdomain 6:1 (8.25..8.5)
                -0x1.e755a2e90fed5p-11, 0x1.0a2ac7c5cca9ap-4, // subdomain 6:2 (8.5..8.75)
                -0x1.ccc263b906921p-11, 0x1.02c3554c1d773p-4, // subdomain 6:3 (8.75..9.0)
                -0x1.b449e1ca6d840p-11, 0x1.f783cd62d82c6p-5, // subdomain 6:4 (9.0..9.25)
                -0x1.9db4fa9a996a4p-11, 0x1.ea3d037a209f6p-5, // subdomain 6:5 (9.25..9.5)
                -0x1.88d363ba5e5f5p-11, 0x1.dda3f567c9b7ap-5, // subdomain 6:6 (9.5..9.75)
                -0x1.757aad4ad1e98p-11, 0x1.d1abb604a6b6fp-5, // subdomain 6:7 (9.75..10.0)
                -0x1.63856e15f97f3p-11, 0x1.c648970cc38a0p-5, // subdomain 6:8 (10.0..10.25)
                -0x1.52d291aaf81b0p-11, 0x1.bb7003c6293edp-5, // subdomain 6:9 (10.25..10.5)
                -0x1.4344c26a378fbp-11, 0x1.b11860c362a18p-5, // subdomain 6:10 (10.5..10.75)
                -0x1.34c1ea9eb4e7ap-11, 0x1.a738eff68002dp-5, // subdomain 6:11 (10.75..11.0)
                -0x1.2732c8b197456p-11, 0x1.9dc9b86c8a3eap-5, // subdomain 6:12 (11.0..11.25)
                -0x1.1a8293434ab7cp-11, 0x1.94c37125d6657p-5, // subdomain 6:13 (11.25..11.5)
                -0x1.0e9eaa8ca4b60p-11, 0x1.8c1f6e96ddb05p-5, // subdomain 6:14 (11.5..11.75)
                -0x1.037654e39e858p-11, 0x1.83d792703353fp-5, // subdomain 6:15 (11.75..12.0)
            };
            // clang-format on
        } // namespace erfcxtab

        namespace derfcxtab
        {
            // clang-format off
            alignas(64) const double coeffs0[144 * 8] = {
                0x1.eeb5403939df2p-78, -0x1.2c64ffcba5108p-70, 0x1.5ae628d0472bap-62, -0x1.7aae51fcf6b7ap-54, 0x1.83acfb47e6136p-46, -0x1.7061d0c15ee88p-38, 0x1.403e24b73286dp-30, -0x1.f27b0fabcb6e6p-23, // subdomain 0:0 (0.5..0.515625)
                0x1.ce68e215dab8cp-78, -0x1.19b54965aa172p-70, 0x1.467170b006897p-62, -0x1.65a5303a11f50p-54, 0x1.6f8c308c1859bp-46, -0x1.5eafa035709b7p-38, 0x1.323396c7502c2p-30, -0x1.dee86dbf4d604p-23, // subdomain 0:1 (0.515625..0.53125)
                0x1.b04fde67da470p-78, -0x1.083d6bd4c1a77p-70, 0x1.3341f64e4fdb3p-62, -0x1.51d9153fc7676p-54, 0x1.5c8922c29848dp-46, -0x1.4de835669367cp-38, 0x1.24d55db701c15p-30, -0x1.cc3105d0b28f3p-23, // subdomain 0:2 (0.53125..0.546875)
                0x1.94424a2c815a5p-78, -0x1.efd05f52760c6p-71, 0x1.21428e18af0d4p-62, -0x1.3f3658616f880p-54, 0x1.4a9306b44df1ap-46, -0x1.3dfe914fea34dp-38, 0x1.181a8f6fe9fdap-30, -0x1.ba4a5a8ffa097p-23, // subdomain 0:3 (0.546875..0.5625)
                0x1.7a1b2934c71dfp-78, -0x1.d143bb8d0fa02p-71, 0x1.105f7cac5cfe7p-62, -0x1.2daa980907342p-54, 0x1.399a1b571ba68p-46, -0x1.2ee678243c515p-38, 0x1.0bfabff76078fp-30, -0x1.a92a79538054ap-23, // subdomain 0:4 (0.5625..0.578125)
                0x1.61b83516dd164p-78, -0x1.b4b046602a36ep-71, 0x1.00865cd2188cfp-62, -0x1.1d24a386c258ap-54, 0x1.298f9887d9963p-46, -0x1.2094653cb471bp-38, 0x1.006dfa0b483d9p-30, -0x1.98c7f271fe4dfp-23, // subdomain 0:5 (0.578125..0.59375)
                0x1.4af9a8992eac3p-78, -0x1.99f3edf40cfc4p-71, 0x1.e34c0ec40e799p-63, -0x1.0d94667171bd9p-54, 0x1.1a659eefde062p-46, -0x1.12fd7fcfccabfp-38, 0x1.ead970673c9e6p-31, -0x1.8919d20f0a74bp-23, // subdomain 0:6 (0.59375..0.609375)
                0x1.35c20f3a4fbf5p-78, -0x1.80ef0465d3459p-71, 0x1.c75cf9d54033ep-63, -0x1.fdd5aaf4ce265p-55, 0x1.0c0f2900f6a51p-46, -0x1.06179062dd1a2p-38, 0x1.d5dfbc8265da1p-31, -0x1.7a1799530e3eap-23, // subdomain 0:7 (0.609375..0.625)
                0x1.21f618810495dp-78, -0x1.6984133dde9d9p-71, 0x1.ad21a1f11065dp-63, -0x1.e233b93375878p-55, 0x1.fcfff9e884a48p-47, -0x1.f3b1edd36d69ep-39, 0x1.c1e1666505cfep-31, -0x1.6bb9380815c6bp-23, // subdomain 0:8 (0.625..0.640625)
                0x1.0f7c6ec694b11p-78, -0x1.5397b240c6cfbp-71, 0x1.947e2dbfe227bp-63, -0x1.c8289cf3c61b7p-55, 0x1.e3593f738f7c9p-47, -0x1.dc7143112c2cfp-39, 0x1.aed1b6f03f9a2p-31, -0x1.5df706954cfbap-23, // subdomain 0:9 (0.640625..0.65625)
                0x1.fc7b22681075bp-79, -0x1.3f1061580df98p-71, 0x1.7d589fd4f7ef4p-63, -0x1.af9baa16de179p-55, 0x1.cb14918fed9d4p-47, -0x1.c65c07db16f54p-39, 0x1.9ca4a7250d17fp-31, -0x1.50c9c05362f49p-23, // subdomain 0:10 (0.65625..0.671875)
                0x1.dc476161a4702p-79, -0x1.2bd665553677fp-71, 0x1.6798b59dc8fb9p-63, -0x1.9875c656ebdc7p-55, 0x1.b41dad76f59b6p-47, -0x1.b1621e6e6f8a0p-39, 0x1.8b4ed6075a815p-31, -0x1.442a7e346ac63p-23, // subdomain 0:11 (0.671875..0.6875)
                0x1.be351f0d994bfp-79, -0x1.19d3a756bfb55p-71, 0x1.5327c8b38f6f9p-63, -0x1.82a14e91cb64bp-55, 0x1.9e618971c854dp-47, -0x1.9d7454bee55e8p-39, 0x1.7ac57f1de4f62p-31, -0x1.3812b1ba29903p-23, // subdomain 0:12 (0.6875..0.703125)
                0x1.a21ec708b448fp-79, -0x1.08f3969b0eef9p-71, 0x1.3ff0b261883a1p-63, -0x1.6e09fdecafd4ap-55, 0x1.89ce410577529p-47, -0x1.8a8456422b354p-39, 0x1.6afe7184e68bcp-31, -0x1.2c7c203613f6ep-23, // subdomain 0:13 (0.703125..0.71875)
                0x1.87e174d409df2p-79, -0x1.f2461922d565ap-72, 0x1.2ddfb135dfc30p-63, -0x1.5a9cd6aee2532p-55, 0x1.7653026f10338p-47, -0x1.78849ea18395bp-39, 0x1.5bf0078a43919p-31, -0x1.2160de4e8a3a7p-23, // subdomain 0:14 (0.71875..0.734375)
                0x1.6f5cc0ea2bb62p-79, -0x1.d4a065f9f5eb3p-72, 0x1.1ce2507680e36p-63, -0x1.48480cb4257a3p-55, 0x1.63dffd5868dfbp-47, -0x1.67686d4211ea1p-39, 0x1.4d911eca8a49bp-31, -0x1.16bb4bc52a02ap-23, // subdomain 0:15 (0.734375..0.75)
                0x1.587291b8d04a9p-79, -0x1.b8d4d8027777cp-72, 0x1.0ce75155f4cfdp-63, -0x1.36faf159a1442p-55, 0x1.526652b02fe46p-47, -0x1.5723b995c9becp-39, 0x1.3fd910c6a9005p-31, -0x1.0c860f7a5019ap-23, // subdomain 0:16 (0.75..0.765625)
                0x1.4306f0306fd9ap-79, -0x1.9ec477e507182p-72, 0x1.fbbd2b8e7f2d1p-64, -0x1.26a5e0ca623d3p-55, 0x1.41d805914ca09p-47, -0x1.47ab2827ce24cp-39, 0x1.32bfabeac7006p-31, -0x1.02bc13aa23098p-23, // subdomain 0:17 (0.765625..0.78125)
                0x1.2effdfa1afb02p-79, -0x1.865269905fcc7p-72, 0x1.df7219a4552bfp-64, -0x1.173a3092739b6p-55, 0x1.3227ed270638bp-47, -0x1.38f40057fa0d7p-39, 0x1.263d2cff2c224p-31, -0x1.f2b104bf744f6p-24, // subdomain 0:18 (0.78125..0.796875)
                0x1.1c4538a64acf8p-79, -0x1.6f63c6217bb31p-72, 0x1.c4d1409782a0fp-64, -0x1.08aa1f6171f38p-55, 0x1.2349a77cb8d84p-47, -0x1.2af422ba20e4bp-39, 0x1.1a4a38fc993a4p-31, -0x1.e0ad843656bb7p-24, // subdomain 0:19 (0.796875..0.8125)
                0x1.0ac086d8808fep-79, -0x1.59df789521babp-72, 0x1.abc047b1d2af5p-64, -0x1.f5d18bce4d465p-56, 0x1.15318d29189f3p-47, -0x1.1da2000e5dee6p-39, 0x1.0edfd73de223fp-31, -0x1.cf64e54efa55ap-24, // subdomain 0:20 (0.8125..0.828125)
                0x1.f4b9d241f4300p-80, -0x1.45ae1d0ce5b3cp-72, 0x1.94268a45ef02ep-64, -0x1.dbd4116c9af9cp-56, 0x1.07d4a5c622979p-47, -0x1.10f490c88215bp-39, 0x1.03f76c09011cfp-31, -0x1.beced43b546d9p-24, // subdomain 0:21 (0.828125..0.84375)
                0x1.d60de8cf2fa1cp-80, -0x1.32b9e275e423cp-72, 0x1.7decfa4fed363p-64, -0x1.c345163632c74p-56, 0x1.f6513a4fcfea1p-48, -0x1.04e34d1d5bffcp-39, 0x1.f31566d4763bbp-32, -0x1.aee3642421123p-24, // subdomain 0:22 (0.84375..0.859375)
                0x1.b95930ef7271ap-80, -0x1.20ee6e6360815p-72, 0x1.68fe052403206p-64, -0x1.ac0f45a8c74c0p-56, 0x1.de47728ab071fp-48, -0x1.f2cc4b1a6c70ep-40, 0x1.df2778b88cc08p-32, -0x1.9f9b09d3a0a48p-24, // subdomain 0:23 (0.859375..0.875)
                0x1.9e7a0dd8b7683p-80, -0x1.1038c2f2faa4ap-72, 0x1.55457a053eb82p-64, -0x1.961e993f5e80bp-56, 0x1.c779a18e2a556p-48, -0x1.dceaf7c6f8472p-40, 0x1.cc19c8716f720p-32, -0x1.90ee96ab8ab13p-24, // subdomain 0:24 (0.875..0.890625)
                0x1.855137f7cf0b2p-80, -0x1.0087269375b1ep-72, 0x1.42b07280f99acp-64, -0x1.816043129d6c4p-56, 0x1.b1d68462625c6p-48, -0x1.c8143948018ccp-40, 0x1.b9e1a4d892951p-32, -0x1.82d733e1d573bp-24, // subdomain 0:25 (0.890625..0.90625)
                0x1.6dc191fd6765ep-80, -0x1.e3921b141b238p-73, 0x1.312d3c6e63218p-64, -0x1.6dc299e6158b3p-56, 0x1.9d4dd73afd650p-48, -0x1.b43a71bdc686dp-40, 0x1.a874e8d5ab72fp-32, -0x1.754e5e0035f2dp-24, // subdomain 0:26 (0.90625..0.921875)
                0x1.57b0011c70be3p-80, -0x1.c7de0a2a3b44dp-73, 0x1.20ab4573eace1p-64, -0x1.5b350689470fbp-56, 0x1.89d045f99ea59p-48, -0x1.a150c0c9096b5p-40, 0x1.97c9f3be4db07p-32, -0x1.684de0a27cccbp-24, // subdomain 0:27 (0.921875..0.9375)
                0x1.4303483b31628p-80, -0x1.add5401ad1413p-73, 0x1.111b07f693ce3p-64, -0x1.49a7f274dc52bp-56, 0x1.774f5dac57db7p-48, -0x1.8f4af8a6e5359p-40, 0x1.87d7a225ec0b6p-32, -0x1.5bcfd2702fdbbp-24, // subdomain 0:28 (0.9375..0.953125)
                0x1.2fa3e5dd4dad9p-80, -0x1.955cc9d356ff3p-73, 0x1.026df9594d3bcp-64, -0x1.390cb78e51753p-56, 0x1.65bd7ef75ac45p-48, -0x1.7e1d93f521972p-40, 0x1.789547176f1bap-32, -0x1.4fce914df91a0p-24, // subdomain 0:29 (0.953125..0.96875)
                0x1.1d7bf491956f8p-80, -0x1.7e5b7a149e488p-73, 0x1.e92cf2e8aa008p-65, -0x1.295590fdc9f7bp-56, 0x1.550dd15a616afp-48, -0x1.6dbdac145d501p-40, 0x1.69faa5b203078p-32, -0x1.4444bec3b7d6dp-24, // subdomain 0:30 (0.96875..0.984375)
                0x1.0c770db295c52p-80, -0x1.68b9ca6a1739cp-73, 0x1.cf0f865ce044fp-65, -0x1.1a758d03475a9p-56, 0x1.4534374353fbap-48, -0x1.5e20f01e1a3bap-40, 0x1.5bffeb2324249p-32, -0x1.392d3c9432a90p-24, // subdomain 0:31 (0.984375..1.0)
                0x1.e9d8b85a4d628p-71, -0x1.4ab0faa50fc23p-64, 0x1.aa99352f55178p-57, -0x1.059e7ad9207dcp-49, 0x1.2ee662084e737p-42, -0x1.480e911dd515fp-35, 0x1.48237ca3af254p-28, -0x1.2955b33169aa2p-21, // subdomain 1:0 (1.0..1.03125)
                0x1.b1acdf4460b9dp-71, -0x1.26a99369f61bcp-64, 0x1.7eb02d2ea8b02p-57, -0x1.d8b4c076bcb94p-50, 0x1.13b12d214a381p-42, -0x1.2cf4fd038fca2p-35, 0x1.2f9599336d978p-28, -0x1.159abd0f50f5dp-21, // subdomain 1:1 (1.03125..1.0625)
                0x1.803bc8086191fp-71, -0x1.06c0d36705aadp-64, 0x1.578c3199d43e5p-57, -0x1.ab5bf11f3c1f0p-50, 0x1.f634da9f98154p-43, -0x1.1448673b45e85p-35, 0x1.190d0aa94a7c7p-28, -0x1.035839d03e410p-21, // subdomain 1:2 (1.0625..1.09375)
                0x1.54ae9dc68c58dp-71, -0x1.d4f171abde2b7p-65, 0x1.34a209b38c12bp-57, -0x1.82a2c4f3f1419p-50, 0x1.c9bb1b4228e51p-43, -0x1.fb9a76fd3956ap-36, 0x1.045baba841474p-28, -0x1.e4de970a56912p-22, // subdomain 1:3 (1.09375..1.125)
                0x1.2e4a57c10026fp-71, -0x1.a2c5dc0619af5p-65, 0x1.1577134e5c420p-57, -0x1.5e09d1897c359p-50, 0x1.a17ac1e65c3dap-43, -0x1.d29bf1ad53556p-36, 0x1.e2af80ff19b0cp-29, -0x1.c587a54b0d5b8p-22, // subdomain 1:4 (1.125..1.15625)
                0x1.0c6c008fb829fp-71, -0x1.763d86f875b06p-65, 0x1.f33e680db4aa3p-58, -0x1.3d20138d2cad7p-50, 0x1.7d06960bf072dp-43, -0x1.ad351453e1d48p-36, 0x1.bfb708a8f3731p-29, -0x1.a87856392301bp-22, // subdomain 1:5 (1.15625..1.1875)
                0x1.dd0b020e6e00bp-72, -0x1.4eaeded29f364p-65, 0x1.c1761c9cbed2dp-58, -0x1.1f813b1533242p-50, 0x1.5bfd06090ccf0p-43, -0x1.8b0f6370d81c3p-36, 0x1.9f89853fb34dbp-29, -0x1.8d81e03f13ec8p-22, // subdomain 1:6 (1.1875..1.21875)
                0x1.a835b8220fd47p-72, -0x1.2b8511095b9afp-65, 0x1.94ece27a74825p-58, -0x1.04d42d8e88275p-50, 0x1.3e06d8c4785adp-43, -0x1.6bdd0a91ce019p-36, 0x1.81e8f736b8c1fp-29, -0x1.74798c258568fp-22, // subdomain 1:7 (1.21875..1.25)
                0x1.797f96e868e03p-72, -0x1.0c3d690a8778ap-65, 0x1.6d0df6c2f9002p-58, -0x1.d9936cde30131p-51, 0x1.22d6069c316a8p-43, -0x1.4f57f222f56b4p-36, 0x1.669d24c9fdb1ap-29, -0x1.5d3853b15653fp-22, // subdomain 1:8 (1.25..1.28125)
                0x1.502b5fc235193p-72, -0x1.e0ca0bcdfd5e5p-66, 0x1.4955d3a13f8d3p-58, -0x1.ae36c1754b0c8p-51, 0x1.0a24b4b01a43fp-43, -0x1.3540ef44d3427p-36, 0x1.4d73078712aa0p-29, -0x1.479a89f07b4cdp-22, // subdomain 1:9 (1.28125..1.3125)
                0x1.2b93911ff115fp-72, -0x1.af2db97f6b757p-66, 0x1.29501b59a36e8p-58, -0x1.8714e8621c757p-51, 0x1.e7689cc746574p-44, -0x1.1d5f0aa097e47p-36, 0x1.363c492ebe42fp-29, -0x1.337f8c3121036p-22, // subdomain 1:10 (1.3125..1.34375)
                0x1.0b2750e37ffb6p-72, -0x1.82f3fb287bea8p-66, 0x1.0c95c630b0c16p-58, -0x1.63be23533ae2fp-51, 0x1.be9972e05d2adp-44, -0x1.077edb92b5967p-36, 0x1.20cece3e80da2p-29, -0x1.20c97abb573fep-22, // subdomain 1:11 (1.34375..1.375)
                0x1.dccf818b80537p-73, -0x1.5b802306f6e5ap-66, 0x1.e59714fbd423cp-59, -0x1.43ced569694aap-51, 0x1.997742ae48951p-44, -0x1.e6e3eaaa94e8fp-37, 0x1.0d044ca31f9acp-29, -0x1.0f5cf87efb9a0p-22, // subdomain 1:12 (1.375..1.40625)
                0x1.a9cb5c25e77b9p-73, -0x1.38481599702cdp-66, 0x1.b740eec9fca23p-59, -0x1.26ee223d1143dp-51, 0x1.77a7aed148e72p-44, -0x1.c21cc83750813p-37, 0x1.f573da827a2c0p-30, -0x1.fe41e1fc4420ep-23, // subdomain 1:13 (1.40625..1.4375)
                0x1.7c7f205049645p-73, -0x1.18d1ffdea2e8bp-66, 0x1.8d997762628dfp-59, -0x1.0cccb6bba3801p-51, 0x1.58d989d92eab9p-44, -0x1.a05c709e17976p-37, 0x1.d39fee4f5470ep-30, -0x1.dffcc7a05b451p-23, // subdomain 1:14 (1.4375..1.46875)
                0x1.543f95c92dbffp-73, -0x1.f964b0c2912f2p-67, 0x1.682126592430cp-59, -0x1.ea476d69eae6ep-52, 0x1.3cc3da926adccp-44, -0x1.815e3ca72c9b1p-37, 0x1.b45306c24dd2ep-30, -0x1.c3c06a56f4817p-23, // subdomain 1:15 (1.46875..1.5)
                0x1.307657915f75cp-73, -0x1.c7144260db53ap-67, 0x1.4666ab180ec55p-59, -0x1.bf679532883e4p-52, 0x1.2324fc8b34ec3p-44, -0x1.64e402d9a8d0ep-37, 0x1.975870f29eda7p-30, -0x1.a96605e651825p-23, // subdomain 1:16 (1.5..1.53125)
                0x1.109f3392e136fp-73, -0x1.9a0ac47920d3cp-67, 0x1.2805455d0519fp-59, -0x1.988892f5a67e9p-52, 0x1.0bc1d983a8c61p-44, -0x1.4ab57106dcd20p-37, 0x1.7c8019c5e2c8ep-30, -0x1.90c9fba912bb1p-23, // subdomain 1:17 (1.53125..1.5625)
                0x1.e88bc3f48a8b1p-74, -0x1.71b3a973b442ep-67, 0x1.0ca350f900fabp-59, -0x1.7544f7ca6b90ep-52, 0x1.ecca71c51e3a3p-45, -0x1.329f776587e8fp-37, 0x1.639e1f0acf2b2p-30, -0x1.79cb8c17fd20ep-23, // subdomain 1:18 (1.5625..1.59375)
                0x1.b60810e91e039p-74, -0x1.4d8b724ce8223p-67, 0x1.e7e1fcc205179p-60, -0x1.5541e91587cecp-52, 0x1.c5be4544b4ab2p-45, -0x1.1c73c34bb3ac9p-37, 0x1.4c8a6b968e19fp-30, -0x1.644c96e85eea1p-23, // subdomain 1:19 (1.59375..1.625)
                0x1.88fefe4d5eb8ap-74, -0x1.2d1da4de8f146p-67, 0x1.bb4e64f00d5f4p-60, -0x1.382df851b354bp-52, 0x1.a208a6f0e70bcp-45, -0x1.080847cf68754p-37, 0x1.37205d3dddebap-30, -0x1.503161184ddd4p-23, // subdomain 1:20 (1.625..1.65625)
                0x1.60d19d35c5897p-74, -0x1.100302d4e102ap-67, 0x1.930d11647ddb4p-60, -0x1.1dc01cae3cd1dp-52, 0x1.815b807b680cfp-45, -0x1.ea6da591e2b6bp-38, 0x1.233e73a371298p-30, -0x1.3d6060535f8e6p-23, // subdomain 1:21 (1.65625..1.6875)
                0x1.3cf3be28e7c8bp-74, -0x1.ebbff20672e04p-68, 0x1.6eace6a46c1b7p-60, -0x1.05b6ca7de3ad7p-52, 0x1.63704db63b218p-45, -0x1.c7b959bbbea73p-38, 0x1.10c606f1ad386p-30, -0x1.2bc20b29fa787p-23, // subdomain 1:22 (1.6875..1.71875)
                0x1.1ce9a5862eb1cp-74, -0x1.bcc67fccb9717p-68, 0x1.4dc8ec1762534p-60, -0x1.dfae49de7ebbbp-53, 0x1.480756654c82ep-45, -0x1.a7b486c65481cp-38, 0x1.ff360b5e8ca44p-31, -0x1.1b40ada32de9cp-23, // subdomain 1:23 (1.71875..1.75)
                0x1.00460a65157bbp-74, -0x1.92894e29970ecp-68, 0x1.3006efa1bdfaap-60, -0x1.b7d88ded8bfe1p-53, 0x1.2ee6fd5465d7fp-45, -0x1.8a25b60eb6f1ep-38, 0x1.df4771e96ea05p-31, -0x1.0bc841bd4f33ep-23, // subdomain 1:24 (1.75..1.78125)
                0x1.cd50ac5af4b88p-75, -0x1.6c880c5f108ddp-68, 0x1.151651ee8f2fcp-60, -0x1.938d41667fd75p-53, 0x1.17db224ee9859p-45, -0x1.6ed894a851d78p-38, 0x1.c1911cb91460ap-31, -0x1.fa8c96f92f923p-24, // subdomain 1:25 (1.78125..1.8125)
                0x1.9f7638f476bc5p-75, -0x1.4a509d365b78ap-68, 0x1.f95decf450f47p-61, -0x1.7276cdd3eca61p-53, 0x1.02b494de0e99ep-45, -0x1.559d7663a411cp-38, 0x1.a5e7ce343ef9dp-31, -0x1.df537061a0170p-24, // subdomain 1:26 (1.8125..1.84375)
                0x1.766587f291ab2p-75, -0x1.2b7d73cc862efp-68, 0x1.cd20a61a323a1p-61, -0x1.5448250390f86p-53, 0x1.de912be202f8cp-46, -0x1.3e48e55cad6acp-38, 0x1.8c23d764c3744p-31, -0x1.c5c58147ae7b5p-24, // subdomain 1:27 (1.84375..1.875)
                0x1.5198b3ec6c8d3p-75, -0x1.0fb4229375bf3p-68, 0x1.a50135c28d3abp-61, -0x1.38bbdc49e12bdp-53, 0x1.bae0cd92d020fp-46, -0x1.28b33cbf4c871p-38, 0x1.7420c7bfa321bp-31, -0x1.adc59bd560095p-24, // subdomain 1:28 (1.875..1.90625)
                0x1.30991627a03c1p-75, -0x1.ed482c9998f41p-69, 0x1.8098044608e8dp-61, -0x1.1f9360d3a02d6p-53, 0x1.9a11c76e894a4p-46, -0x1.14b84d94faa37p-38, 0x1.5dbd246f96a75p-31, -0x1.9738c13379238p-24, // subdomain 1:29 (1.90625..1.9375)
                0x1.12fd78e669c40p-75, -0x1.c00aef216c14fp-69, 0x1.5f882dd1841f6p-61, -0x1.089642294cd3fp-53, 0x1.7be452b3c923ep-46, -0x1.02370c9040352p-38, 0x1.48da2663f8739p-31, -0x1.8205f4463a35fp-24, // subdomain 1:30 (1.9375..1.96875)
                0x1.f0d1049e2d43bp-76, -0x1.97305c47f1814p-69, 0x1.417e5a2e76526p-61, -0x1.e72320be05280p-54, 0x1.601e86824aa59p-46, -0x1.e2228fced0c56p-39, 0x1.355b7e79baab6p-31, -0x1.6e16105786c41p-24, // subdomain 1:31 (1.96875..2.0)
                0x1.ab197c9d3ae46p-66, -0x1.613c19a9ce302p-60, 0x1.196cc9faaf686p-53, -0x1.aeb1d761305b5p-47, 0x1.3a866b649015dp-40, -0x1.b350d9974bf8ep-34, 0x1.1a813829f1801p-27, -0x1.525d3407e77a8p-21, // subdomain 2:0 (2.0..2.0625)
                0x1.5dbeed0287efcp-66, -0x1.24c0ec5e26c95p-60, 0x1.d856707861ef3p-54, -0x1.6e285cf80d1fcp-47, 0x1.0f0ce02795f35p-40, -0x1.7c8742c4ec3b1p-34, 0x1.f560b6c76927ap-28, -0x1.311583a188f19p-21, // subdomain 2:1 (2.0625..2.125)
                0x1.1f0f92982d975p-66, -0x1.e657b2161a13ap-61, 0x1.8d3d92f4c8f9bp-54, -0x1.37f1275e9de10p-47, 0x1.d41ac01613eb3p-41, -0x1.4d4481b812becp-34, 0x1.bdb3e86133cf9p-28, -0x1.1386b63e56770p-21, // subdomain 2:2 (2.125..2.1875)
                0x1.d849cce855fa8p-67, -0x1.94dcab8d9fdfap-61, 0x1.4ecb79fc95d62p-54, -0x1.0a4c8f4fe6d48p-47, 0x1.950069dfc3846p-41, -0x1.246b3124f0220p-34, 0x1.8ce5e45c5edefp-28, -0x1.f274d33e0cb26p-22, // subdomain 2:3 (2.1875..2.25)
                0x1.8561e0b21c4e9p-67, -0x1.51c2f420ffa36p-61, 0x1.1ac181b67890fp-54, -0x1.c795917b4bb81p-48, 0x1.5f14349457813p-41, -0x1.010b10567c6dfp-34, 0x1.62097f29f23ecp-28, -0x1.c394990185420p-22, // subdomain 2:4 (2.25..2.3125)
                0x1.41bb9c86acfaep-67, -0x1.1a61991698263p-61, 0x1.de97b9f4342b6p-55, -0x1.8679ec2e60b09p-48, 0x1.30e8b35ee14f4p-41, -0x1.c4b2e3777bbdcp-35, 0x1.3c5411037d768p-28, -0x1.99bb88d124e9ep-22, // subdomain 2:5 (2.3125..2.375)
                0x1.0a68c00439316p-67, -0x1.d9263ab0567e7p-62, 0x1.95d9b61d51650p-55, -0x1.4f52491b715c5p-48, 0x1.094cbe000f76fp-41, -0x1.8f55ff3c632fep-35, 0x1.1b181791da24cp-28, -0x1.744fe34b7b14bp-22, // subdomain 2:6 (2.375..2.4375)
                0x1.ba21e1afaa1b9p-68, -0x1.8d362a0cfb922p-62, 0x1.58d8b79bbed05p-55, -0x1.208134a49d573p-48, 0x1.ce81ed1521f87p-42, -0x1.60dd2f57b860cp-35, 0x1.fb8174a3b5e49p-29, -0x1.52cc03d390311p-22, // subdomain 2:7 (2.4375..2.5)
                0x1.6fa559a9ff603p-68, -0x1.4e227ea62a7b5p-62, 0x1.259524b2cfdc3p-55, -0x1.f15ffe0fb2d78p-49, 0x1.93de70523857ap-42, -0x1.385214777acc4p-35, 0x1.c79c14904806ep-29, -0x1.34bb8920836b7p-22, // subdomain 2:8 (2.5..2.5625)
                0x1.3255c66a3416cp-68, -0x1.19a1dff6d978ep-62, 0x1.f4d59bcb7b152p-56, -0x1.ad840d68428d9p-49, 0x1.614723a4fc372p-42, -0x1.14e327f4a76d1p-35, 0x1.99a3b98d1c4e9p-29, -0x1.19b8ebd7f72b4p-22, // subdomain 2:9 (2.5625..2.625)
                0x1.ff86f238c7189p-69, -0x1.dbaf73c3672dbp-63, 0x1.abffb87d8b921p-56, -0x1.7394ca57054f3p-49, 0x1.358d04837625fp-42, -0x1.ebbbc2026b1ffp-36, 0x1.70da4dafe40ddp-29, -0x1.016b7185bebf4p-22, // subdomain 2:10 (2.625..2.6875)
                0x1.abedbae7e3af2p-69, -0x1.927e6235bdb61p-63, 0x1.6e6e2221370e5p-56, -0x1.42070e5dcb5a3p-49, 0x1.0fb0725f99a80p-42, -0x1.b553aa8544f26p-36, 0x1.4c9b874111eecp-29, -0x1.d70ada84e1e5ep-23, // subdomain 2:11 (2.6875..2.75)
                0x1.66b0ee1b9a0f1p-69, -0x1.5534a32e42e0dp-63, 0x1.3a49579f54f8fp-56, -0x1.1790a239d64a2p-49, 0x1.ddb2b5e081429p-43, -0x1.8589511ce79c6p-36, 0x1.2c5927b0d47abp-29, -0x1.af8583a259cddp-23, // subdomain 2:12 (2.75..2.8125)
                0x1.2d3badd30d570p-69, -0x1.21c92c2117153p-63, 0x1.0e0adcd4018dcp-56, -0x1.e63a0e1e8f309p-50, 0x1.a4a1848639f56p-43, -0x1.5b7d4f4cf4d1fp-36, 0x1.0f97cf4e062e7p-29, -0x1.8bcf3608caeafp-23, // subdomain 2:13 (2.8125..2.875)
                0x1.fae9e45a849c3p-70, -0x1.ed205fe7cb1b3p-64, 0x1.d0dd525bb833dp-57, -0x1.a7888c9955460p-50, 0x1.72f6cee47edd6p-43, -0x1.366fe0d4b2336p-36, 0x1.ebd89ad34d136p-30, -0x1.6b7eabcccb32fp-23, // subdomain 2:14 (2.875..2.9375)
                0x1.ab4eda56fb099p-70, -0x1.a453735621b90p-64, 0x1.90cea2a518397p-57, -0x1.7186cd2d6f3dcp-50, 0x1.47aaa59502a00p-43, -0x1.15bc1e95676e9p-36, 0x1.bdf2b1811cbb7p-30, -0x1.4e37210a03232p-23, // subdomain 2:15 (2.9375..3.0)
                0x1.68dc96faa25c5p-70, -0x1.66e74cf685928p-64, 0x1.5a28c74cc65b4p-57, -0x1.42ec6e521c0c5p-50, 0x1.21dc47934f507p-43, -0x1.f1a7f87f2f70ap-37, 0x1.94db3ee6ee0fcp-30, -0x1.33a6b6cbe672dp-23, // subdomain 2:16 (3.0..3.0625)
                0x1.314b9d912e276p-70, -0x1.32fca398689b9p-64, 0x1.2b75161976d87p-57, -0x1.1aa464ab57e6fp-50, 0x1.00cbfa2bb237cp-43, -0x1.be79d02e5b467p-37, 0x1.7004414e7c257p-30, -0x1.1b850f3b28656p-23, // subdomain 2:17 (3.0625..3.125)
                0x1.02bdeb582c7fbp-70, -0x1.0706d37b4e2b6p-64, 0x1.037a41e84daf8p-57, -0x1.ef890fe805feap-51, 0x1.c7abcb3282c7fp-44, -0x1.9119e3bac3e11p-37, 0x1.4ef188734ee2ep-30, -0x1.05921a99bc086p-23, // subdomain 2:18 (3.125..3.1875)
                0x1.b75644bddfd53p-71, -0x1.c379fc349b61fp-65, 0x1.c263a3ffd617ep-58, -0x1.b30d2a01ff356p-51, 0x1.94db924ac17cfp-44, -0x1.68d0c5406150cp-37, 0x1.31364dc276a89p-30, -0x1.e32a1bd4115e4p-24, // subdomain 2:19 (3.1875..3.25)
                0x1.75a124f3df33ep-71, -0x1.841be615feed1p-65, 0x1.877ef9563efe5p-58, -0x1.7e84c36175e22p-51, 0x1.6836c38e3bc94p-44, -0x1.44ff660bf6c60p-37, 0x1.167325486b00fp-30, -0x1.beb6fa720a5eep-24, // subdomain 2:20 (3.25..3.3125)
                0x1.3e48b11f6bc5fp-71, -0x1.4e2caee24a066p-65, 0x1.54d415fceafa1p-58, -0x1.50d2185efeb46p-51, 0x1.40eec9c62dddap-44, -0x1.251ba1a386418p-37, 0x1.fca87388d236ap-31, -0x1.9d7128df94c07p-24, // subdomain 2:21 (3.3125..3.375)
                0x1.0f95ec3676957p-71, -0x1.20313ab93480ep-65, 0x1.292afa33e0989p-58, -0x1.2901f238321c3p-51, 0x1.1e51ec6469be5p-44, -0x1.08ad4de4b1e93p-37, 0x1.d11f910baf53bp-31, -0x1.7f08d4af88e71p-24, // subdomain 2:22 (3.375..3.4375)
                0x1.d03b8da874e67p-72, -0x1.f1d9c9d6959d2p-66, 0x1.037c6852fedb9p-58, -0x1.0644f88f0eb97p-51, 0x1.ff8e10ee69571p-45, -0x1.de97763eeb75bp-38, 0x1.a9c9a70ecce9ap-31, -0x1.6336d446eef3ep-24, // subdomain 2:23 (3.4375..3.5)
                0x1.8d66e34fef82cp-72, -0x1.aead7b8147c88p-66, 0x1.c5d3ca4cdd4d5p-59, -0x1.cfd439a561cb1p-52, 0x1.c993e99721c70p-45, -0x1.b137252c6c46dp-38, 0x1.8634037acd82cp-31, -0x1.49bba0b4486b8p-24, // subdomain 2:24 (3.5..3.5625)
                0x1.54ba7dc3df299p-72, -0x1.75216cc21504ap-66, 0x1.8d6e2eee8fe88p-59, -0x1.9ab3de003b01dp-52, 0x1.99d0ea4dc0f65p-45, -0x1.889a0964f8557p-38, 0x1.65f94d9d76b30p-31, -0x1.325e70f5cc7bcp-24, // subdomain 2:25 (3.5625..3.625)
                0x1.249727770399ap-72, -0x1.43c0bd9a3a824p-66, 0x1.5c88900440443p-59, -0x1.6c2567c22df51p-52, 0x1.6f7f4aa4ecbfap-45, -0x1.643475484c27ep-38, 0x1.48bfd8f9964eap-31, -0x1.1cec721f8c57fp-24, // subdomain 2:26 (3.625..3.6875)
                0x1.f746127466363p-73, -0x1.1951d7382800dp-66, 0x1.32139d83c1773p-59, -0x1.434a380cb2f92p-52, 0x1.49f39151ae1eap-45, -0x1.438c023c94dc9p-38, 0x1.2e3831c051e33p-31, -0x1.0938187b6935fp-24, // subdomain 2:27 (3.6875..3.75)
                0x1.b17b2a320d8b9p-73, -0x1.e9991dac4bd16p-67, 0x1.0d2887b9f6dd5p-59, -0x1.1f630104e4a8fp-52, 0x1.2898dc2c117c4p-45, -0x1.26354718750cap-38, 0x1.161bdab92fddfp-31, -0x1.ee310c97cf6fdp-25, // subdomain 2:28 (3.75..3.8125)
                0x1.75eb42558d996p-73, -0x1.aaa3fd4ee6715p-67, 0x1.da05487e879c3p-60, -0x1.ff963c8bc26a9p-53, 0x1.0aedb656e2ac0p-45, -0x1.0bd1e0f2651d3p-38, 0x1.002c359baf92ep-31, -0x1.ccd20aa8c42ddp-25, // subdomain 2:29 (3.8125..3.875)
                0x1.4302928a5f1e3p-73, -0x1.744b8babd0275p-67, 0x1.a1f43c665d932p-60, -0x1.c7e94f09cef63p-53, 0x1.e102cc7c81e61p-46, -0x1.e81d842b25674p-39, 0x1.d8631fd0ff48fp-32, -0x1.ae112160ac414p-25, // subdomain 2:30 (3.875..3.9375)
                0x1.176e85c3ba5a7p-73, -0x1.4550deac86afcp-67, 0x1.70ff07e65f34dp-60, -0x1.96ca3de22b3f8p-53, 0x1.b1e3418ed421ap-46, -0x1.bd45798d0ecd6p-39, 0x1.b3f49e634f01cp-32, -0x1.91b2d9b6dbe28p-25, // subdomain 2:31 (3.9375..4.0)
                0x1.c305bf984ee7cp-64, -0x1.0a7eca06be6b2p-58, 0x1.32b98b1a4e77ep-52, -0x1.579e45174c476p-46, 0x1.747d768cc2c02p-40, -0x1.84b1bf80ba517p-34, 0x1.832cb4f7bfc90p-28, -0x1.6b2a8af969b1bp-22, // subdomain 3:0 (4.0..4.125)
                0x1.54585cbb6c673p-64, -0x1.9a317af41f849p-59, 0x1.e1c1f3b376107p-53, -0x1.137f2482f1e43p-46, 0x1.31130202bf2b3p-40, -0x1.4567bfae29166p-34, 0x1.4b8ccff972aa9p-28, -0x1.3e594c89e72cep-22, // subdomain 3:1 (4.125..4.25)
                0x1.0232f534e431ap-64, -0x1.3d4af850f1303p-59, 0x1.7c27b0ed7477ep-53, -0x1.bbc07991851b0p-47, 0x1.f5cf99cbbe91bp-41, -0x1.11771b6735650p-34, 0x1.1ce53a7673b32p-28, -0x1.17e7cc3fe1cc3p-22, // subdomain 3:2 (4.25..4.375)
                0x1.89c989401097cp-65, -0x1.ed47e5fd9caffp-60, 0x1.2d5fdd882ffecp-53, -0x1.66f06d1051d3fp-47, 0x1.9e60191872b4ap-41, -0x1.cd53e1f16e9aep-35, 0x1.eb3be53449646p-29, -0x1.eda2d7a968bdap-23, // subdomain 3:3 (4.375..4.5)
                0x1.2dcb6fd519c3cp-65, -0x1.81448f7442ed8p-60, 0x1.dffcb87d59720p-54, -0x1.238eb8b9a04fcp-47, 0x1.5781c0ecb24edp-41, -0x1.86819d21a89bap-35, 0x1.a8db794ca805ep-29, -0x1.b47edf2093a6bp-23, // subdomain 3:4 (4.5..4.625)
                0x1.d0d58478ac9c6p-66, -0x1.2e4b73f5eee6ep-60, 0x1.7fe3de19711dbp-54, -0x1.db932ee9015c4p-48, 0x1.1dd41d410b78bp-41, -0x1.4bb09223e6e71p-35, 0x1.7093932171a91p-29, -0x1.83006c33dba0dp-23, // subdomain 3:5 (4.625..4.75)
                0x1.67a9f0394e222p-66, -0x1.dc80a8de0c570p-61, 0x1.3452f7ebcca1cp-54, -0x1.85642f2689477p-48, 0x1.dd65470529924p-42, -0x1.1aa99c68bcadfp-35, 0x1.40b1f6e94cd88p-29, -0x1.58007e00db54ap-23, // subdomain 3:6 (4.75..4.875)
                0x1.178f96784dbf2p-66, -0x1.792d9e9661a79p-61, 0x1.f1483ef47f0fbp-55, -0x1.40096e2c47e23p-48, 0x1.9012e04746780p-42, -0x1.e34c97cab2772p-36, 0x1.17d3de2a3a9bdp-29, -0x1.328934dffaf4ap-23, // subdomain 3:7 (4.875..5.0)
                0x1.b48601af60e5ep-67, -0x1.2bcf795d3378dp-61, 0x1.929a705b3a60dp-55, -0x1.080056805755ep-48, 0x1.5069012e409ddp-42, -0x1.9e71f8ad5cb48p-36, 0x1.e9ac9c50b661dp-30, -0x1.11cce846e7874p-23, // subdomain 3:8 (5.0..5.125)
                0x1.5645ff88a9cb6p-67, -0x1.de906e7401e3ap-62, 0x1.4731f5dcc617dp-55, -0x1.b51a2db5951d8p-49, 0x1.1bcd7fc74d59cp-42, -0x1.647499a859a3ap-36, 0x1.ad93749474f99p-30, -0x1.ea3df5c24823ap-24, // subdomain 3:9 (5.125..5.25)
                0x1.0d7e52800eb8dp-67, -0x1.7f7444464d700p-62, 0x1.0ae50e42e4273p-55, -0x1.6b188ec602e55p-49, 0x1.e05bc2cc69f68p-43, -0x1.3375632647961p-36, 0x1.79d06c329b53fp-30, -0x1.b7dc1f5a90e3dp-24, // subdomain 3:10 (5.25..5.375)
                0x1.aa1923d139ba6p-68, -0x1.346c006592bc0p-62, 0x1.b4f9c63159015p-56, -0x1.2e9fcc85035c5p-49, 0x1.97c40ab103d18p-43, -0x1.09ee5587fb6d2p-36, 0x1.4d1b4da7a3881p-30, -0x1.8b7eaf241d3d2p-24, // subdomain 3:11 (5.375..5.5)
                0x1.522dd5c51597bp-68, -0x1.f1fc114707a3cp-63, 0x1.66f6d0de94e76p-56, -0x1.fa1315e292de0p-50, 0x1.5b2b00e416496p-43, -0x1.cd419dd164427p-37, 0x1.26630bca57a91p-30, -0x1.64557d616686ep-24, // subdomain 3:12 (5.5..5.625)
                0x1.0d6d25280ec61p-68, -0x1.9379fc3a4a79ep-63, 0x1.27df8d6fd4575p-56, -0x1.a8792528324c6p-50, 0x1.286c1e923da1dp-43, -0x1.910f2025dbf38p-37, 0x1.04c3bdb8fe145p-30, -0x1.41af99d199f3dp-24, // subdomain 3:13 (5.625..5.75)
                0x1.aee5642414c4cp-69, -0x1.480c17ccc2000p-63, 0x1.e95613f79a199p-57, -0x1.651b96512f430p-50, 0x1.fb97ae95620ecp-44, -0x1.5d97552f7896cp-37, 0x1.cefd1ab8b781fp-31, -0x1.22f62cdb23823p-24, // subdomain 3:14 (5.75..5.875)
                0x1.59cf8e64e5662p-69, -0x1.0b9ff40b80220p-63, 0x1.95eec063a14bbp-57, -0x1.2d50c130d8309p-50, 0x1.b3c4762eeb819p-44, -0x1.31776852053a1p-37, 0x1.9be662e253656p-31, -0x1.07a8413b82ec1p-24, // subdomain 3:15 (5.875..6.0)
                0x1.167fd30059504p-69, -0x1.b61a4accf4df8p-64, 0x1.51c7e115a7eabp-57, -0x1.fded5bb5f97b9p-51, 0x1.771462d0a9464p-44, -0x1.0b89826fc9dbep-37, 0x1.6f35451593637p-31, -0x1.deae9411bdb9dp-25, // subdomain 3:16 (6.0..6.125)
                0x1.c21bbd2a3ff83p-70, -0x1.67bc817f1ed1ep-64, 0x1.19e91ccb64fabp-57, -0x1.b0ad064b6a924p-51, 0x1.43a8851c83103p-44, -0x1.d5b3163e09682p-38, 0x1.480597b294f7ap-31, -0x1.b34883841403ep-25, // subdomain 3:17 (6.125..6.25)
                0x1.6cef9756484d7p-70, -0x1.284eab9c831d4p-64, 0x1.d7ed9f5bddfb8p-58, -0x1.701c92dc0e511p-51, 0x1.17f849309165bp-44, -0x1.9d383e0637ae1p-38, 0x1.25960bce1323cp-31, -0x1.8c7a8c615cdfdp-25, // subdomain 3:18 (6.25..6.375)
                0x1.28d5754ea5a49p-70, -0x1.e999db4a408c9p-65, 0x1.8c201447e2b13p-58, -0x1.39ff1322fe719p-51, 0x1.e581bf8123804p-45, -0x1.6c4e8bc86a147p-38, 0x1.074267db6b5cbp-31, -0x1.69b734da3b7c1p-25, // subdomain 3:19 (6.375..6.5)
                0x1.e465180947a29p-71, -0x1.95af1641c9c3dp-65, 0x1.4d689738e19ddp-58, -0x1.0c83564c4d2d3p-51, 0x1.a5f00df7c7144p-45, -0x1.41d91997490b9p-38, 0x1.d8fd9209d4464p-32, -0x1.4a83f96fe8ba5p-25, // subdomain 3:20 (6.5..6.625)
                0x1.8c70e7ccd3f93p-71, -0x1.511c9b44f39c1p-65, 0x1.195e462b0829ep-58, -0x1.cc5ce476a3308p-52, 0x1.6f82f46724d53p-45, -0x1.1ce8786ab28dep-38, 0x1.a9a772d43ae81p-32, -0x1.2e7682b7ccdf2p-25, // subdomain 3:21 (6.625..6.75)
                0x1.456bb807a39d7p-71, -0x1.18e8e4b87065ep-65, 0x1.dc214d7292943p-59, -0x1.8b959cbf3d02bp-52, 0x1.40ccf18a75d80p-45, -0x1.f965c6d9b3fa3p-39, 0x1.7fb5e6f4582ddp-32, -0x1.15324d26dd7fdp-25, // subdomain 3:22 (6.75..6.875)
                0x1.0be546ae978efp-71, -0x1.d56c338b337efp-66, 0x1.93ddbdc3b0195p-59, -0x1.54b564102c336p-52, 0x1.189de9401b290p-45, -0x1.c11bc1779a57dp-39, 0x1.5a796502b552fp-32, -0x1.fccd5ea4d63f9p-26, // subdomain 3:23 (6.875..7.0)
                0x1.ba51f8144dbf0p-72, -0x1.89414c1cceb72p-66, 0x1.57688fd2b7bf0p-59, -0x1.261b253af211cp-52, 0x1.ebf044aafb8ecp-46, -0x1.8fd2116a9c080p-39, 0x1.395aed7a0b42dp-32, -0x1.d39a5d530db6dp-26, // subdomain 3:24 (7.0..7.125)
                0x1.6e2826634b362p-72, -0x1.4a4aa0f764340p-66, 0x1.24b24be7d90fcp-59, -0x1.fcde2a8c97f92p-53, 0x1.b00eb8ac446aep-46, -0x1.649376cb9f400p-39, 0x1.1bd850794376dp-32, -0x1.ae502b1b0ba82p-26, // subdomain 3:25 (7.125..7.25)
                0x1.2feaf75c5dfe1p-72, -0x1.161a20be8d3bep-66, 0x1.f41a5ee48bb1fp-60, -0x1.b92a112b54a4cp-53, 0x1.7c33239c30d84p-46, -0x1.3e8f2e3080391p-39, 0x1.01810f8fcb414p-32, -0x1.8c82800bbb7cdp-26, // subdomain 3:26 (7.25..7.375)
                0x1.f9d46cc643c23p-73, -0x1.d57534b4f1cb7p-67, 0x1.ac3422ea79b3bp-60, -0x1.7f42f197c3bf4p-53, 0x1.4f322436a2e8dp-46, -0x1.1d1317fcc0f3fp-39, 0x1.d3e780a47e69fp-33, -0x1.6dd22c608a99cp-26, // subdomain 3:27 (7.375..7.5)
                0x1.a602b2f09f6a5p-73, -0x1.8d2e1a96881e8p-67, 0x1.6f7268713a1e6p-60, -0x1.4da159eb2d769p-53, 0x1.280f5f2a74e8fp-46, -0x1.ff0dc0083618ap-40, 0x1.a9b7b3069aacdp-33, -0x1.51eb5a3073e81p-26, // subdomain 3:28 (7.5..7.625)
                0x1.60f3a8362b79fp-73, -0x1.50ce82398e0ddp-67, 0x1.3bfc57fdf0a3cp-60, -0x1.22ff28bba6e48p-53, 0x1.05f59948d2774p-46, -0x1.cacfedc247809p-40, 0x1.83dfb15a7d164p-33, -0x1.38840fc53c5f7p-26, // subdomain 3:29 (7.625..7.75)
                0x1.27e7fcd0a26f9p-73, -0x1.1e40ed3e5d705p-67, 0x1.104c661118124p-60, -0x1.fc97f13bf41edp-54, 0x1.d06074f40531fp-47, -0x1.9c8c6ead3136cp-40, 0x1.61df629ffa335p-33, -0x1.215ae88265547p-26, // subdomain 3:30 (7.75..7.875)
                0x1.f1543045c6b3dp-74, -0x1.e7a52dc661818p-68, 0x1.d641a175b5e64p-61, -0x1.bd476e3905105p-54, 0x1.9c4beae2f2dd3p-47, -0x1.7382326eb5061p-40, 0x1.4346f2dd12d1fp-33, -0x1.0c35fbf2f0126p-26, // subdomain 3:31 (7.875..8.0)
                0x1.80f737563700cp-64, -0x1.811634e8a72b1p-59, 0x1.7a8da28234763p-53, -0x1.6e031baafe23fp-47, 0x1.59fad40042fdbp-41, -0x1.3e58345c5e88ep-35, 0x1.1af273226b4b1p-29, -0x1.dfaf7e5dd0656p-24, // subdomain 4:0 (8.0..8.25)
                0x1.1385f45f46c37p-64, -0x1.1ae5743037d3fp-59, 0x1.1d8ce5bb51adbp-53, -0x1.1b8bbef8f173bp-47, 0x1.135b975f75371p-41, -0x1.0460288b51ad3p-35, 0x1.dbd043b72865dp-30, -0x1.9ec46e722efdep-24, // subdomain 4:1 (8.25..8.5)
                0x1.8dbac7f87e788p-65, -0x1.a2ed66f3f822cp-60, 0x1.b1ee3b804f7dbp-54, -0x1.ba3faab438ec0p-48, 0x1.b8f1850306704p-42, -0x1.ac35976ece0a7p-36, 0x1.91f1951d78417p-30, -0x1.680ed6f1bdba1p-24, // subdomain 4:2 (8.5..8.75)
                0x1.2163b2db00beap-65, -0x1.3884615922b58p-60, 0x1.4bfe8599dfca1p-54, -0x1.5b13dc9214b64p-48, 0x1.631331eac0b79p-42, -0x1.61e851540aae5p-36, 0x1.550cdb09fb1b1p-30, -0x1.39bdde89ee176p-24, // subdomain 4:3 (8.75..9.0)
                0x1.a85ec0dd3b863p-66, -0x1.d59de5085db0ap-61, 0x1.ff5e667460e69p-55, -0x1.1208a17201966p-48, 0x1.1f7d067610201p-42, -0x1.25eaa1f0f03dbp-36, 0x1.229a92c6b77abp-30, -0x1.125c6164e495ep-24, // subdomain 4:4 (9.0..9.25)
                0x1.3972d6f0fed0fp-66, -0x1.63430b2f954b2p-61, 0x1.8c50b9a725601p-55, -0x1.b338189ce3304p-49, 0x1.d3f42386b4626p-43, -0x1.ea7190072fa36p-37, 0x1.f13a6fef53912p-31, -0x1.e17a5570a5508p-25, // subdomain 4:5 (9.25..9.5)
                0x1.d24d980a8ada6p-67, -0x1.0e834ec9f6022p-61, 0x1.34ffffafad1c4p-55, -0x1.5b822c4a92bdap-49, 0x1.7ebd8c7116b20p-43, -0x1.9afd90f8f98cap-37, 0x1.ab02f6c1fdfaep-31, -0x1.a7d793a289c82p-25, // subdomain 4:6 (9.5..9.75)
                0x1.5d304f77dc1c3p-67, -0x1.9e8ada5fe0df6p-62, 0x1.e49f0a6bcf9e3p-56, -0x1.16eeee9dd3e81p-49, 0x1.3a86eb8a69a18p-43, -0x1.59db165573a7ap-37, 0x1.700ca21da0568p-31, -0x1.7641b8c98eb49p-25, // subdomain 4:7 (9.75..10.0)
                0x1.072d81acd34b4p-67, -0x1.3f8881ad70cd7p-62, 0x1.7e1ff926b97edp-56, -0x1.c20840fd43027p-50, 0x1.03a4145e30777p-43, -0x1.24362f9b547dap-37, 0x1.3e558ef63043dp-31, -0x1.4b71ebd624dfdp-25, // subdomain 4:8 (10.0..10.25)
                0x1.8f2993c90d45bp-68, -0x1.ef6d9ad368c1fp-63, 0x1.2ee5d49c3e420p-56, -0x1.6cca02de82e9bp-50, 0x1.ae858697d5cd5p-44, -0x1.efab9ace9255bp-38, 0x1.1440563f577e1p-31, -0x1.265b07e690435p-25, // subdomain 4:9 (10.25..10.5)
                0x1.3080f36d7fda3p-68, -0x1.82309077cd925p-63, 0x1.e29ff6a519795p-57, -0x1.290e30c479917p-50, 0x1.6669ad9e0f4efp-44, -0x1.a5f0643c6bf6bp-38, 0x1.e0fd85009f1e1p-32, -0x1.061e703949fabp-25, // subdomain 4:10 (10.5..10.75)
                0x1.d33c224adaf45p-69, -0x1.2ea0e2615d70fp-63, 0x1.825d30ebbb652p-57, -0x1.e5ef58b5bdcb4p-51, 0x1.2b9126c82f961p-44, -0x1.6870294d0dae1p-38, 0x1.a402ba430ff17p-32, -0x1.d40678a68faf7p-26, // subdomain 4:11 (10.75..11.0)
                0x1.686e05aadc6b1p-69, -0x1.dcb3b9e3a17acp-64, 0x1.36be771e459bfp-57, -0x1.8f240f29b1d8fp-51, 0x1.f6ac58ea40b39p-45, -0x1.34f06a2f3d4bfp-38, 0x1.6fd54e1792024p-32, -0x1.a2de6a321dd4dp-26, // subdomain 4:12 (11.0..11.25)
                0x1.1780de2809512p-69, -0x1.79485587ff6d9p-64, 0x1.f616173dc8a61p-58, -0x1.492f106940bedp-51, 0x1.a748dd27762e6p-45, -0x1.09a780742b95cp-38, 0x1.430a3a801408ep-32, -0x1.77c2782888ab8p-26, // subdomain 4:13 (11.25..11.5)
                0x1.b3b0b45c4f8f0p-70, -0x1.2bff641c2dbfcp-64, 0x1.975e9139f4250p-58, -0x1.108cc9f95046ap-51, 0x1.65af894d5bf83p-45, -0x1.ca493b5ed732dp-39, 0x1.1c76c49c402cap-32, -0x1.51d99cafef551p-26, // subdomain 4:14 (11.5..11.75)
                0x1.553ba08a6b9e9p-70, -0x1.df3d9a14d7e2cp-65, 0x1.4be365b977d82p-58, -0x1.c5042cbc7b3a3p-52, 0x1.2f45ea36a44c4p-45, -0x1.8c7a7a26fa927p-39, 0x1.f6483b8dae608p-33, -0x1.306ccd44bdd45p-26, // subdomain 4:15 (11.75..12.0)
            };

            alignas(64) const double coeffs1[144 * 4] = {
                0x1.4f1b8ac55ce2fp-15, -0x1.6ab030260bd40p-8, 0x1.03a4819c57725p-1, 0, // subdomain 0:0 (0.5..0.515625)
                0x1.43a7d0a04911cp-15, -0x1.60657102402fap-8, 0x1.fc1ce9c776aabp-2, 0, // subdomain 0:1 (0.515625..0.53125)
                0x1.38a6ef363169fp-15, -0x1.567480dfd309ap-8, 0x1.f141bcaed92d2p-2, 0, // subdomain 0:2 (0.53125..0.546875)
                0x1.2e13e2ea62412p-15, -0x1.4cd9dd315404dp-8, 0x1.e6b4bb9bbedbcp-2, 0, // subdomain 0:3 (0.546875..0.5625)
                0x1.23e9e56c240d4p-15, -0x1.43922a8edddf5p-8, 0x1.dc7341b19e96ap-2, 0, // subdomain 0:4 (0.5625..0.578125)
                0x1.1a246a8e024f3p-15, -0x1.3a9a32d8780dep-8, 0x1.d27ac458d8a75p-2, 0, // subdomain 0:5 (0.578125..0.59375)
                0x1.10bf1d499b579p-15, -0x1.31eee3710a317p-8, 0x1.c8c8d21b9b795p-2, 0, // subdomain 0:6 (0.59375..0.609375)
                0x1.07b5dced5f3e2p-15, -0x1.298d4b9087a4ep-8, 0x1.bf5b119092f5ep-2, 0, // subdomain 0:7 (0.609375..0.625)
                0x1.fe0974e39ab26p-16, -0x1.21729aac0d84fp-8, 0x1.b62f4052ade75p-2, 0, // subdomain 0:8 (0.625..0.640625)
                0x1.ed4febe7cba49p-16, -0x1.199c1ef2c0821p-8, 0x1.ad4332054ea6ap-2, 0, // subdomain 0:9 (0.640625..0.65625)
                0x1.dd37f8a5519b6p-16, -0x1.120743de498abp-8, 0x1.a494cf6445a3cp-2, 0, // subdomain 0:10 (0.65625..0.671875)
                0x1.cdbac9db3c4f9p-16, -0x1.0ab190d5e10d2p-8, 0x1.9c22155efd1ebp-2, 0, // subdomain 0:11 (0.671875..0.6875)
                0x1.bed1df050065cp-16, -0x1.0398a7e2e8350p-8, 0x1.93e9143e45a7cp-2, 0, // subdomain 0:12 (0.6875..0.703125)
                0x1.b0770450828e7p-16, -0x1.f97488ec3c826p-9, 0x1.8be7eed43ac99p-2, 0, // subdomain 0:13 (0.703125..0.71875)
                0x1.a2a44ecb7f008p-16, -0x1.ec2874791bbebp-9, 0x1.841cd9b5be948p-2, 0, // subdomain 0:14 (0.71875..0.734375)
                0x1.955418c5279a3p-16, -0x1.df48e806b6f04p-9, 0x1.7c861a7d11b04p-2, 0, // subdomain 0:15 (0.734375..0.75)
                0x1.8880fe6103473p-16, -0x1.d2d1e5535316fp-9, 0x1.75220715141f5p-2, 0, // subdomain 0:16 (0.75..0.765625)
                0x1.7c25da584693ap-16, -0x1.c6bf983fa76dap-9, 0x1.6def050cc1093p-2, 0, // subdomain 0:17 (0.765625..0.78125)
                0x1.703dc2e707c48p-16, -0x1.bb0e54e697931p-9, 0x1.66eb88f27daf5p-2, 0, // subdomain 0:18 (0.78125..0.796875)
                0x1.64c406e2d8981p-16, -0x1.afba95ccdce4ap-9, 0x1.601615b6d92c0p-2, 0, // subdomain 0:19 (0.796875..0.8125)
                0x1.59b42af8751d2p-16, -0x1.a4c0fa275d36cp-9, 0x1.596d3c165fca0p-2, 0, // subdomain 0:20 (0.8125..0.828125)
                0x1.4f09e70e5a0d5p-16, -0x1.9a1e4436ff10ep-9, 0x1.52ef9a0a29a2cp-2, 0, // subdomain 0:21 (0.828125..0.84375)
                0x1.44c123ca35e63p-16, -0x1.8fcf57b8dc77cp-9, 0x1.4c9bda3ed0d21p-2, 0, // subdomain 0:22 (0.84375..0.859375)
                0x1.3ad5f83746c5fp-16, -0x1.85d13869c53cbp-9, 0x1.4670b3917fe85p-2, 0, // subdomain 0:23 (0.859375..0.875)
                0x1.3144a78bd3f52p-16, -0x1.7c21089c10c23p-9, 0x1.406ce892cd605p-2, 0, // subdomain 0:24 (0.875..0.890625)
                0x1.28099f0c0d226p-16, -0x1.72bc07decd427p-9, 0x1.3a8f470f1cc60p-2, 0, // subdomain 0:25 (0.890625..0.90625)
                0x1.1f217408b2c94p-16, -0x1.699f91b567ec1p-9, 0x1.34d6a79c41e38p-2, 0, // subdomain 0:26 (0.90625..0.921875)
                0x1.1688e1f80324bp-16, -0x1.60c91c5ef5ae7p-9, 0x1.2f41ed2c25c67p-2, 0, // subdomain 0:27 (0.921875..0.9375)
                0x1.0e3cc8a77d7bap-16, -0x1.583637ac5057cp-9, 0x1.29d004a430b74p-2, 0, // subdomain 0:28 (0.9375..0.953125)
                0x1.063a2a8522ae2p-16, -0x1.4fe48be446c49p-9, 0x1.247fe4793f578p-2, 0, // subdomain 0:29 (0.953125..0.96875)
                0x1.fcfc55fddb6bdp-17, -0x1.47d1d8b529588p-9, 0x1.1f508c4fec0b3p-2, 0, // subdomain 0:30 (0.96875..0.984375)
                0x1.ee0c19eea2ce8p-17, -0x1.3ffbf43305d77p-9, 0x1.1a4104a0f8998p-2, 0, // subdomain 0:31 (0.984375..1.0)
                0x1.d896e4f966ea7p-15, -0x1.34a899c7ce572p-8, 0x1.12e35755d082ep-2, 0, // subdomain 1:0 (1.0..1.03125)
                0x1.bda6b15798d28p-15, -0x1.2650df14e479bp-8, 0x1.0978011ad5cc5p-2, 0, // subdomain 1:1 (1.03125..1.0625)
                0x1.a47e1ebe5a41fp-15, -0x1.18c96ff487903p-8, 0x1.007c1e0967a43p-2, 0, // subdomain 1:2 (1.0625..1.09375)
                0x1.8cfb5a93fe055p-15, -0x1.0c049896a4e8fp-8, 0x1.efd2c6911df50p-3, 0, // subdomain 1:3 (1.09375..1.125)
                0x1.76ff559343e40p-15, -0x1.ffeb510d824ddp-9, 0x1.df73df0356426p-3, 0, // subdomain 1:4 (1.125..1.15625)
                0x1.626d852377373p-15, -0x1.e921bb2046a3fp-9, 0x1.cfd08636a4ffbp-3, 0, // subdomain 1:5 (1.15625..1.1875)
                0x1.4f2baaa458c31p-15, -0x1.d396a0fe067d6p-9, 0x1.c0de722afc48fp-3, 0, // subdomain 1:6 (1.1875..1.21875)
                0x1.3d21a013e3919p-15, -0x1.bf35c72320a44p-9, 0x1.b29400f19fee3p-3, 0, // subdomain 1:7 (1.21875..1.25)
                0x1.2c39298469a26p-15, -0x1.abec64e0644afp-9, 0x1.a4e82c983448ap-3, 0, // subdomain 1:8 (1.25..1.28125)
                0x1.1c5dcae81afccp-15, -0x1.99a90710ea50fp-9, 0x1.97d28009056eep-3, 0, // subdomain 1:9 (1.28125..1.3125)
                0x1.0d7ca1c3ff4fcp-15, -0x1.885b75536fe5bp-9, 0x1.8b4b0ccb67f66p-3, 0, // subdomain 1:10 (1.3125..1.34375)
                0x1.ff0884cfe0a18p-16, -0x1.77f4998bdd692p-9, 0x1.7f4a619108a5dp-3, 0, // subdomain 1:11 (1.34375..1.375)
                0x1.e4c930a6c8266p-16, -0x1.6866697777e97p-9, 0x1.73c9817ec369fp-3, 0, // subdomain 1:12 (1.375..1.40625)
                0x1.cc1d92eff02bdp-16, -0x1.59a3d22374325p-9, 0x1.68c1dc2130a6cp-3, 0, // subdomain 1:13 (1.40625..1.4375)
                0x1.b4ea3825a16d7p-16, -0x1.4ba0a51a51024p-9, 0x1.5e2d45fe87370p-3, 0, // subdomain 1:14 (1.4375..1.46875)
                0x1.9f15bdba444f6p-16, -0x1.3e5187209197ap-9, 0x1.5405f1b8bd226p-3, 0, // subdomain 1:15 (1.46875..1.5)
                0x1.8a88a6b562eb1p-16, -0x1.31abe05d37af5p-9, 0x1.4a4669b3fcedap-3, 0, // subdomain 1:16 (1.5..1.53125)
                0x1.772d3425c86ddp-16, -0x1.25a5cdcdc0076p-9, 0x1.40e98a3694a13p-3, 0, // subdomain 1:17 (1.53125..1.5625)
                0x1.64ef410c4001ap-16, -0x1.1a3613e870f3cp-9, 0x1.37ea7bf87920cp-3, 0, // subdomain 1:18 (1.5625..1.59375)
                0x1.53bc216bb2960p-16, -0x1.0f54125288ddap-9, 0x1.2f44af1956c7ep-3, 0, // subdomain 1:19 (1.59375..1.625)
                0x1.43828433a723ap-16, -0x1.04f7b89252dedp-9, 0x1.26f3d674f1c3ap-3, 0, // subdomain 1:20 (1.625..1.65625)
                0x1.343257c18fdb9p-16, -0x1.f632f74ec2c94p-10, 0x1.1ef3e34e4f610p-3, 0, // subdomain 1:21 (1.65625..1.6875)
                0x1.25bcb0bbf06ffp-16, -0x1.e36498ea60e60p-10, 0x1.1741014ac75dfp-3, 0, // subdomain 1:22 (1.6875..1.71875)
                0x1.1813b3105bbdfp-16, -0x1.d1771ddc81088p-10, 0x1.0fd792b6b2e23p-3, 0, // subdomain 1:23 (1.71875..1.75)
                0x1.0b2a7ce29cb14p-16, -0x1.c05e23debab54p-10, 0x1.08b42d0df7965p-3, 0, // subdomain 1:24 (1.75..1.78125)
                0x1.fdea268050367p-17, -0x1.b00e0f787befdp-10, 0x1.01d395c329b39p-3, 0, // subdomain 1:25 (1.78125..1.8125)
                0x1.e6d0a0de84e14p-17, -0x1.a07bfe2000ff0p-10, 0x1.f6657e80e2216p-4, 0, // subdomain 1:26 (1.8125..1.84375)
                0x1.d0f3a76a44780p-17, -0x1.919db96a7363ep-10, 0x1.e99d8c3b8864ep-4, 0, // subdomain 1:27 (1.84375..1.875)
                0x1.bc3fe4e3a4c65p-17, -0x1.8369ab34ce84ap-10, 0x1.dd49dd10d46b0p-4, 0, // subdomain 1:28 (1.875..1.90625)
                0x1.a8a354e7714fap-17, -0x1.75d6d2b122dfep-10, 0x1.d16543b013d70p-4, 0, // subdomain 1:29 (1.90625..1.9375)
                0x1.960d2aceecd06p-17, -0x1.68dcba459b652p-10, 0x1.c5ead89bb4f4ep-4, 0, // subdomain 1:30 (1.9375..1.96875)
                0x1.846dba96c7cb7p-17, -0x1.5c736e2c48839p-10, 0x1.bad5f5f5a1e56p-4, 0, // subdomain 1:31 (1.96875..2.0)
                0x1.6bad73b1c05f8p-15, -0x1.4ad4be02541f5p-9, 0x1.aaeb6de6cd87ap-4, 0, // subdomain 2:0 (2.0..2.0625)
                0x1.4d8c069f38a39p-15, -0x1.350f1adcd7194p-9, 0x1.96eed1c6ff575p-4, 0, // subdomain 2:1 (2.0625..2.125)
                0x1.325ba8eea1f54p-15, -0x1.21138efb1a061p-9, 0x1.84400062608bdp-4, 0, // subdomain 2:2 (2.125..2.1875)
                0x1.19c8b4b2a2573p-15, -0x1.0eb5b56134518p-9, 0x1.72c3c2647ea03p-4, 0, // subdomain 2:3 (2.1875..2.25)
                0x1.0389d74a04249p-15, -0x1.fb9c1dab9cf2cp-10, 0x1.62617ebe608b9p-4, 0, // subdomain 2:4 (2.25..2.3125)
                0x1.debd59a120dc4p-16, -0x1.dc72d03a2ae68p-10, 0x1.5302f136393eap-4, 0, // subdomain 2:5 (2.3125..2.375)
                0x1.ba1d1eefce2b9p-16, -0x1.bfb0a9a336cd7p-10, 0x1.4493e9ea28bd0p-4, 0, // subdomain 2:6 (2.375..2.4375)
                0x1.98cf25698ed91p-16, -0x1.a51d77be6c443p-10, 0x1.37021494b8ad4p-4, 0, // subdomain 2:7 (2.4375..2.5)
                0x1.7a7b4704ea762p-16, -0x1.8c86e63446ae8p-10, 0x1.2a3cc68db2093p-4, 0, // subdomain 2:8 (2.5..2.5625)
                0x1.5ed3902f84e4bp-16, -0x1.75bfcf9f135b3p-10, 0x1.1e34d2a71729cp-4, 0, // subdomain 2:9 (2.5625..2.625)
                0x1.4592f4a533e89p-16, -0x1.609fa4fe60506p-10, 0x1.12dc622592727p-4, 0, // subdomain 2:10 (2.625..2.6875)
                0x1.2e7c31c0f0a3dp-16, -0x1.4d01e86776513p-10, 0x1.0826d22e707fap-4, 0, // subdomain 2:11 (2.6875..2.75)
                0x1.1958d79651c62p-16, -0x1.3ac5b85251b41p-10, 0x1.fc112a380e0ffp-5, 0, // subdomain 2:12 (2.75..2.8125)
                0x1.05f8732e51711p-16, -0x1.29cd6945013b3p-10, 0x1.e8ee2e79b21c7p-5, 0, // subdomain 2:13 (2.8125..2.875)
                0x1.e85faa43d00d4p-17, -0x1.19fe2bf2722cdp-10, 0x1.d6d14d24c7eacp-5, 0, // subdomain 2:14 (2.875..2.9375)
                0x1.c7b0e1105b009p-17, -0x1.0b3fbe2802c0ap-10, 0x1.c5a8ba656f761p-5, 0, // subdomain 2:15 (2.9375..3.0)
                0x1.a99f9d9513b1fp-17, -0x1.faf84a419ad88p-11, 0x1.b5641bf859ce5p-5, 0, // subdomain 2:16 (3.0..3.0625)
                0x1.8dee2aea593b0p-17, -0x1.e13ee010761e2p-11, 0x1.a5f466aa0682bp-5, 0, // subdomain 2:17 (3.0625..3.125)
                0x1.74653711946b7p-17, -0x1.c92f0340edd32p-11, 0x1.974bbf736e360p-5, 0, // subdomain 2:18 (3.125..3.1875)
                0x1.5cd3185ff2cb1p-17, -0x1.b2a7c04c1eab0p-11, 0x1.895d5fc93f070p-5, 0, // subdomain 2:19 (3.1875..3.25)
                0x1.470b2a12bd7f4p-17, -0x1.9d8b1552c08f0p-11, 0x1.7c1d7cc039a54p-5, 0, // subdomain 2:20 (3.25..3.3125)
                0x1.32e53cf500837p-17, -0x1.89bda6625b2a9p-11, 0x1.6f8130b4d8460p-5, 0, // subdomain 2:21 (3.3125..3.375)
                0x1.203d196f34696p-17, -0x1.77267a2189f48p-11, 0x1.637e672e76076p-5, 0, // subdomain 2:22 (3.375..3.4375)
                0x1.0ef210b985a80p-17, -0x1.65aebddf6068cp-11, 0x1.580bcabee64fep-5, 0, // subdomain 2:23 (3.4375..3.5)
                0x1.fdcd3675c27d6p-18, -0x1.55419023d7849p-11, 0x1.4d20b4a7006aep-5, 0, // subdomain 2:24 (3.5..3.5625)
                0x1.e00004c8c6e0dp-18, -0x1.45cbd0fc50ae2p-11, 0x1.42b51e0f3cea8p-5, 0, // subdomain 2:25 (3.5625..3.625)
                0x1.c44c29280da61p-18, -0x1.373bf75859dc3p-11, 0x1.38c192a942fdap-5, 0, // subdomain 2:26 (3.625..3.6875)
                0x1.aa85c378291c1p-18, -0x1.2981eae0914c2p-11, 0x1.2f3f249449a68p-5, 0, // subdomain 2:27 (3.6875..3.75)
                0x1.9284ff705a935p-18, -0x1.1c8ee1c44b612p-11, 0x1.2627616294eb1p-5, 0, // subdomain 2:28 (3.75..3.8125)
                0x1.7c25aabbd2557p-18, -0x1.1055420adf019p-11, 0x1.1d7448222f70ep-5, 0, // subdomain 2:29 (3.8125..3.875)
                0x1.6746d6ede61e9p-18, -0x1.04c88603973ffp-11, 0x1.1520404e60166p-5, 0, // subdomain 2:30 (3.875..3.9375)
                0x1.53ca85dbe20d5p-18, -0x1.f3ba46f713ddfp-12, 0x1.0d26119250fbfp-5, 0, // subdomain 2:31 (3.9375..4.0)
                0x1.38ed98d33685cp-16, -0x1.d5271052870b4p-11, 0x1.01ccafd71bfc1p-5, 0, // subdomain 3:0 (4.0..4.125)
                0x1.19095097aa09ep-16, -0x1.b012d4ba44e8ap-11, 0x1.e774e0c5f564cp-6, 0, // subdomain 3:1 (4.125..4.25)
                0x1.fa13d7458608bp-17, -0x1.8ebb3c41f48fbp-11, 0x1.cd831a6b740b8p-6, 0, // subdomain 3:2 (4.25..4.375)
                0x1.c8c61aecdbee3p-17, -0x1.70acb1646e895p-11, 0x1.b58bf67a385f3p-6, 0, // subdomain 3:3 (4.375..4.5)
                0x1.9d3d1c400ec12p-17, -0x1.5583bb6d69c7fp-11, 0x1.9f5e138ccd942p-6, 0, // subdomain 3:4 (4.5..4.625)
                0x1.76b134e4b7388p-17, -0x1.3cea7815ecde2p-11, 0x1.8acdd807fd07ep-6, 0, // subdomain 3:5 (4.625..4.75)
                0x1.547920d8a5872p-17, -0x1.269684e327cf0p-11, 0x1.77b4a9eed63d8p-6, 0, // subdomain 3:6 (4.75..4.875)
                0x1.3604e2c751251p-17, -0x1.1247434205207p-11, 0x1.65f04531f3666p-6, 0, // subdomain 3:7 (4.875..5.0)
                0x1.1ad997cf1691fp-17, -0x1.ff88cd3653073p-12, 0x1.55622b5951f3bp-6, 0, // subdomain 3:8 (5.0..5.125)
                0x1.028e0b7bf1c34p-17, -0x1.ddb97de2768e7p-12, 0x1.45ef285897fc9p-6, 0, // subdomain 3:9 (5.125..5.25)
                0x1.d98fcd5ff531dp-18, -0x1.beca6a77f4d25p-12, 0x1.377ee91ea83d4p-6, 0, // subdomain 3:10 (5.25..5.375)
                0x1.b272b7d789629p-18, -0x1.a26fe18d96552p-12, 0x1.29fba10cf878fp-6, 0, // subdomain 3:11 (5.375..5.5)
                0x1.8f3e72df8e582p-18, -0x1.88673d00025f7p-12, 0x1.1d51bc00fe192p-6, 0, // subdomain 3:12 (5.5..5.625)
                0x1.6f7e9cd1dd266p-18, -0x1.7075a8eeb340fp-12, 0x1.116f9b0042596p-6, 0, // subdomain 3:13 (5.625..5.75)
                0x1.52ce28237c6c9p-18, -0x1.5a6719b012c15p-12, 0x1.064559eb0983ep-6, 0, // subdomain 3:14 (5.75..5.875)
                0x1.38d51d41dcc2fp-18, -0x1.460d69078ffa5p-12, 0x1.f78939bb5783cp-7, 0, // subdomain 3:15 (5.875..6.0)
                0x1.2146b86518523p-18, -0x1.333f9446b3d63p-12, 0x1.e3c0c8417a15cp-7, 0, // subdomain 3:16 (6.0..6.125)
                0x1.0bdfd3667ca86p-18, -0x1.21d9161dd12dcp-12, 0x1.d119cb6c875d5p-7, 0, // subdomain 3:17 (6.125..6.25)
                0x1.f0cb1d428ed31p-19, -0x1.11b957c7635e0p-12, 0x1.bf7ed770bd46fp-7, 0, // subdomain 3:18 (6.25..6.375)
                0x1.cd485c87fbf33p-19, -0x1.02c335f51a5a5p-12, 0x1.aedc6dbab6490p-7, 0, // subdomain 3:19 (6.375..6.5)
                0x1.acdc46c529625p-19, -0x1.e9b92cfea3cf2p-13, 0x1.9f20c91f0bab0p-7, 0, // subdomain 3:20 (6.5..6.625)
                0x1.8f36734a6ab53p-19, -0x1.cfdc18abc3bd3p-13, 0x1.903bb039e8f5dp-7, 0, // subdomain 3:21 (6.625..6.75)
                0x1.740fc760191c5p-19, -0x1.b7c50f3235057p-13, 0x1.821e4d2bc1ef9p-7, 0, // subdomain 3:22 (6.75..6.875)
                0x1.5b2942ca55e21p-19, -0x1.a14e202f19ff6p-13, 0x1.74bb09ff2b3a8p-7, 0, // subdomain 3:23 (6.875..7.0)
                0x1.444af7ec3623dp-19, -0x1.8c551168a2d28p-13, 0x1.6805711dc16efp-7, 0, // subdomain 3:24 (7.0..7.125)
                0x1.2f4328d1688b1p-19, -0x1.78baf4e90843dp-13, 0x1.5bf2114fe8b15p-7, 0, // subdomain 3:25 (7.125..7.25)
                0x1.1be58382c894bp-19, -0x1.6663cc47c2512p-13, 0x1.507664d5ade2dp-7, 0, // subdomain 3:26 (7.25..7.375)
                0x1.0a0a78f52b9f5p-19, -0x1.55363752ed125p-13, 0x1.4588bb36c1c40p-7, 0, // subdomain 3:27 (7.375..7.5)
                0x1.f31d554703a1fp-20, -0x1.451b2c8f99e95p-13, 0x1.3b202574db8d1p-7, 0, // subdomain 3:28 (7.5..7.625)
                0x1.d4a4d70b89e92p-20, -0x1.35fdba418e3c5p-13, 0x1.313464572cd5ap-7, 0, // subdomain 3:29 (7.625..7.75)
                0x1.b872a336822e8p-20, -0x1.27cacedb8feb9p-13, 0x1.27bdd88f4494dp-7, 0, // subdomain 3:30 (7.75..7.875)
                0x1.9e53ab5aec31cp-20, -0x1.1a7107e26666dp-13, 0x1.1eb574800d4adp-7, 0, // subdomain 3:31 (7.875..8.0)
                0x1.7aa60aa713a76p-18, -0x1.07df2a8541e56p-12, 0x1.11e93f59f33b9p-7, 0, // subdomain 4:0 (8.0..8.25)
                0x1.50cf1c4df279fp-18, -0x1.e3173af5b3ec8p-13, 0x1.0221654dd3f7bp-7, 0, // subdomain 4:1 (8.25..8.5)
                0x1.2c8fa94c29b01p-18, -0x1.bb4efa1509dacp-13, 0x1.e755a2e90fed5p-8, 0, // subdomain 4:2 (8.5..8.75)
                0x1.0d09474aff514p-18, -0x1.97c0fdbde918dp-13, 0x1.ccc263b906921p-8, 0, // subdomain 4:3 (8.75..9.0)
                0x1.e3094d83f02f2p-19, -0x1.77e1f5fcd0731p-13, 0x1.b449e1ca6d840p-8, 0, // subdomain 4:4 (9.0..9.25)
                0x1.b2d3e3b79d9d7p-19, -0x1.5b3b733d35b82p-13, 0x1.9db4fa9a996a4p-8, 0, // subdomain 4:5 (9.25..9.5)
                0x1.88759800368c2p-19, -0x1.41685ae40b3ccp-13, 0x1.88d363ba5e5f5p-8, 0, // subdomain 4:6 (9.5..9.75)
                0x1.631b23df4b129p-19, -0x1.2a120712e6accp-13, 0x1.757aad4ad1e98p-8, 0, // subdomain 4:7 (9.75..10.0)
                0x1.421324360a48cp-19, -0x1.14edee3cf140dp-13, 0x1.63856e15f97f3p-8, 0, // subdomain 4:8 (10.0..10.25)
                0x1.24c80cb4be7bdp-19, -0x1.01bbb71bdef16p-13, 0x1.52d291aaf81b0p-8, 0, // subdomain 4:9 (10.25..10.5)
                0x1.0abb4a5812a48p-19, -0x1.e08747286b531p-14, 0x1.4344c26a378fbp-8, 0, // subdomain 4:10 (10.5..10.75)
                0x1.e702a95f62fffp-20, -0x1.c0aa83610bf6ap-14, 0x1.34c1ea9eb4e7ap-8, 0, // subdomain 4:11 (10.75..11.0)
                0x1.bd7cf896d2a39p-20, -0x1.a38caac9f1144p-14, 0x1.2732c8b197456p-8, 0, // subdomain 4:12 (11.0..11.25)
                0x1.98489f823a4bbp-20, -0x1.88e3e114140fbp-14, 0x1.1a8293434ab7cp-8, 0, // subdomain 4:13 (11.25..11.5)
                0x1.76def25f6766dp-20, -0x1.706f613f2be8ep-14, 0x1.0e9eaa8ca4b60p-8, 0, // subdomain 4:14 (11.5..11.75)
                0x1.58cbf938306d1p-20, -0x1.59f6372a02e34p-14, 0x1.037654e39e858p-8, 0, // subdomain 4:15 (11.75..12.0)
            };
            // clang-format on
        } // namespace derfcxtab

        namespace anchor
        {
            // clang-format off
            const Anchor anchors[n_zeros] = {
                { // zero 1, |z_n| = 2.4086, tau = 0.0830351, alpha = 2.0928, eta = 11.33 u
                    0x1.fdd0c5610a35ep+0, 0x1.fdd0c5610a35fp+0,
                    -0x1.5ad4d62887d7ep+0, -0x1.5ad4d62887d7dp+0,
                    0x1.c3dc016db42d4p-8,
                    {
                        { // anchor xmym
                            0x1.9a3168b5c2604p-53, -0x1.af3fb10c47151p-55,
                            -0x1.4f690cf365872p-51, 0x1.20dd750429b70p+0,
                            -0x1.875b79336387cp+0, -0x1.1fa1f22822670p+1,
                            0x1.03caad5921654p+2, 0x1.b3638ce7fdc07p-1,
                            -0x1.ed42ca3b34ebfp+1, 0x1.8366e49317e34p+1,
                            -0x1.8d9c74ceb121bp-3, -0x1.35b95a40ece9fp+2,
                            0x1.cca95d9657614p+1, 0x1.0ed8c78a47859p+1,
                            -0x1.67da84ac108f3p+1, 0x1.9262b08693b97p+0,
                            -0x1.09a634402364ep-5, -0x1.21c31f23c5646p+1,
                            0x1.5215e5d8143e6p+0, 0x1.491927a9ee81bp-1,
                            -0x1.632abeadae4e1p-1, 0x1.1bf2e38c431fbp-1,
                            -0x1.012348e3eb384p-3, -0x1.f446703794da5p-2,
                            0x1.120681e1c050bp-2, 0x1.52ee0610ed75fp-5,
                            -0x1.24069622137ccp-4, 0x1.e4677206570aep-4,
                            -0x1.4e84f34ac3052p-5, -0x1.b51467b296ee3p-5,
                            0x1.eb4e35808ef5fp-6, -0x1.262e3b3c0ab56p-7,
                        },
                        { // anchor xmyp
                            -0x1.4f1302a5221b1p-54, -0x1.af3fb10c47154p-55,
                            0x1.dfb5aeb9be8cep-52, 0x1.20dd750429b6dp+0,
                            -0x1.875b79336387fp+0, -0x1.1fa1f2282266ap+1,
                            0x1.03caad5921651p+2, 0x1.b3638ce7fdbe8p-1,
                            -0x1.ed42ca3b34eb3p+1, 0x1.8366e49317e33p+1,
                            -0x1.8d9c74ceb1280p-3, -0x1.35b95a40ece99p+2,
                            0x1.cca95d965760ep+1, 0x1.0ed8c78a47850p+1,
                            -0x1.67da84ac108eap+1, 0x1.9262b08693b96p+0,
                            -0x1.09a6344023707p-5, -0x1.21c31f23c5640p+1,
                            0x1.5215e5d8143e1p+0, 0x1.491927a9ee80dp-1,
                            -0x1.632abeadae4d6p-1, 0x1.1bf2e38c431f8p-1,
                            -0x1.012348e3eb388p-3, -0x1.f446703794d98p-2,
                            0x1.120681e1c0505p-2, 0x1.52ee0610ed741p-5,
                            -0x1.24069622137c0p-4, 0x1.e4677206570a4p-4,
                            -0x1.4e84f34ac304ep-5, -0x1.b51467b296ed5p-5,
                            0x1.eb4e35808ef51p-6, -0x1.262e3b3c0ab58p-7,
                        },
                        { // anchor xpym
                            0x1.9a3168b5c25fcp-53, 0x1.d5eafdc541a87p-53,
                            -0x1.6b62431364874p-50, 0x1.20dd750429b6cp+0,
                            -0x1.875b793363870p+0, -0x1.1fa1f2282266ep+1,
                            0x1.03caad5921650p+2, 0x1.b3638ce7fdc1fp-1,
                            -0x1.ed42ca3b34ebfp+1, 0x1.8366e49317e27p+1,
                            -0x1.8d9c74ceb116ep-3, -0x1.35b95a40ece9cp+2,
                            0x1.cca95d965760ap+1, 0x1.0ed8c78a4785fp+1,
                            -0x1.67da84ac108f3p+1, 0x1.9262b08693b85p+0,
                            -0x1.09a63440234d1p-5, -0x1.21c31f23c5643p+1,
                            0x1.5215e5d8143dfp+0, 0x1.491927a9ee826p-1,
                            -0x1.632abeadae4e3p-1, 0x1.1bf2e38c431f0p-1,
                            -0x1.012348e3eb36bp-3, -0x1.f446703794da3p-2,
                            0x1.120681e1c0507p-2, 0x1.52ee0610ed790p-5,
                            -0x1.24069622137d5p-4, 0x1.e4677206570a2p-4,
                            -0x1.4e84f34ac3044p-5, -0x1.b51467b296ee8p-5,
                            0x1.eb4e35808ef5ep-6, -0x1.262e3b3c0ab39p-7,
                        },
                        { // anchor xpyp
                            -0x1.4f1302a5221b0p-54, 0x1.d5eafdc541a81p-53,
                            -0x1.2f0143ad08823p-52, 0x1.20dd750429b69p+0,
                            -0x1.875b793363872p+0, -0x1.1fa1f22822668p+1,
                            0x1.03caad592164dp+2, 0x1.b3638ce7fdc01p-1,
                            -0x1.ed42ca3b34eb3p+1, 0x1.8366e49317e27p+1,
                            -0x1.8d9c74ceb11d4p-3, -0x1.35b95a40ece96p+2,
                            0x1.cca95d9657604p+1, 0x1.0ed8c78a47855p+1,
                            -0x1.67da84ac108eap+1, 0x1.9262b08693b84p+0,
                            -0x1.09a634402358bp-5, -0x1.21c31f23c563dp+1,
                            0x1.5215e5d8143dap+0, 0x1.491927a9ee818p-1,
                            -0x1.632abeadae4d9p-1, 0x1.1bf2e38c431edp-1,
                            -0x1.012348e3eb36fp-3, -0x1.f446703794d96p-2,
                            0x1.120681e1c0501p-2, 0x1.52ee0610ed773p-5,
                            -0x1.24069622137c9p-4, 0x1.e467720657098p-4,
                            -0x1.4e84f34ac303fp-5, -0x1.b51467b296ed9p-5,
                            0x1.eb4e35808ef50p-6, -0x1.262e3b3c0ab3bp-7,
                        },
                    },
                },
                { // zero 2, |z_n| = 3.4615, tau = 0.0577788, alpha = 2.1067, eta = 11.38 u
                    0x1.5877923be93f2p+1, 0x1.5877923be93f3p+1,
                    -0x1.16a96850ca058p+1, -0x1.16a96850ca057p+1,
                    0x1.b591edfcd6689p-9,
                    {
                        { // anchor xmym
                            0x1.7e1459820dbcfp-52, -0x1.05d0a791fb08ep-51,
                            0x1.beef067d7fb66p-53, 0x1.20dd750429b7fp+0,
                            -0x1.3a6fa4f9e5c09p+1, -0x1.84b08130e3f59p+1,
                            0x1.1a1097772a221p+3, 0x1.216573d694464p+0,
                            -0x1.7b9c8b74e697dp+3, 0x1.32f1cca12a75ep+3,
                            0x1.c84891698f4edp-1, -0x1.51ba54f74b9fbp+4,
                            0x1.278f8277fa464p+4, 0x1.06260d072b41ep+4,
                            -0x1.8a63e5d6eac62p+4, 0x1.3b149e9fc5711p+2,
                            0x1.2928780568d49p+3, -0x1.4d2f1b585d3fcp+4,
                            0x1.3ff471d541264p+3, 0x1.fb433c2f18565p+3,
                            -0x1.c481efecb20f2p+3, -0x1.c195165119b30p-7,
                            0x1.46d2f7cbc3e3bp+2, -0x1.0f2110ced4d6cp+3,
                            0x1.92004ce031fcbp+1, 0x1.69f2a2d321437p+2,
                            -0x1.fd740ba459208p+1, 0x1.c81649339be3ep-7,
                            0x1.13c0baf0dd5a9p+0, -0x1.068ab59cb3d3dp+1,
                            0x1.7aa6eeddcc02fp-1, 0x1.0bfa978bc7798p+0,
                        },
                        { // anchor xmyp
                            -0x1.874d210c8b648p-53, -0x1.05d0a791fb087p-51,
                            0x1.a09f7198bbf02p-49, 0x1.20dd750429b75p+0,
                            -0x1.3a6fa4f9e5c0dp+1, -0x1.84b08130e3f3ep+1,
                            0x1.1a1097772a218p+3, 0x1.216573d694405p+0,
                            -0x1.7b9c8b74e6962p+3, 0x1.32f1cca12a760p+3,
                            0x1.c84891698f364p-1, -0x1.51ba54f74b9edp+4,
                            0x1.278f8277fa460p+4, 0x1.06260d072b409p+4,
                            -0x1.8a63e5d6eac4dp+4, 0x1.3b149e9fc5736p+2,
                            0x1.2928780568d25p+3, -0x1.4d2f1b585d3f1p+4,
                            0x1.3ff471d541264p+3, 0x1.fb433c2f18542p+3,
                            -0x1.c481efecb20dbp+3, -0x1.c195165116304p-7,
                            0x1.46d2f7cbc3e19p+2, -0x1.0f2110ced4d63p+3,
                            0x1.92004ce031fcap+1, 0x1.69f2a2d32141dp+2,
                            -0x1.fd740ba4591ebp+1, 0x1.c81649339cd53p-7,
                            0x1.13c0baf0dd589p+0, -0x1.068ab59cb3d32p+1,
                            0x1.7aa6eeddcc029p-1, 0x1.0bfa978bc7783p+0,
                        },
                        { // anchor xpym
                            0x1.7e1459820dbc7p-52, 0x1.b0ccd722eaea5p-55,
                            -0x1.1e80b4920dc46p-49, 0x1.20dd750429b72p+0,
                            -0x1.3a6fa4f9e5befp+1, -0x1.84b08130e3f55p+1,
                            0x1.1a1097772a215p+3, 0x1.216573d6944b0p+0,
                            -0x1.7b9c8b74e697cp+3, 0x1.32f1cca12a744p+3,
                            0x1.c84891698f6a8p-1, -0x1.51ba54f74b9efp+4,
                            0x1.278f8277fa44fp+4, 0x1.06260d072b422p+4,
                            -0x1.8a63e5d6eac59p+4, 0x1.3b149e9fc56bep+2,
                            0x1.2928780568d5fp+3, -0x1.4d2f1b585d3eap+4,
                            0x1.3ff471d541241p+3, 0x1.fb433c2f18565p+3,
                            -0x1.c481efecb20e4p+3, -0x1.c19516511f863p-7,
                            0x1.46d2f7cbc3e4ep+2, -0x1.0f2110ced4d5bp+3,
                            0x1.92004ce031f97p+1, 0x1.69f2a2d321437p+2,
                            -0x1.fd740ba4591f9p+1, 0x1.c81649339a187p-7,
                            0x1.13c0baf0dd5bfp+0, -0x1.068ab59cb3d2dp+1,
                            0x1.7aa6eeddcc005p-1, 0x1.0bfa978bc779cp+0,
                        },
                        { // anchor xpyp
                            -0x1.874d210c8b628p-53, 0x1.b0ccd722eae77p-55,
                            0x1.98bf327b58c09p-51, 0x1.20dd750429b69p+0,
                            -0x1.3a6fa4f9e5bf2p+1, -0x1.84b08130e3f3bp+1,
                            0x1.1a1097772a20cp+3, 0x1.216573d694451p+0,
                            -0x1.7b9c8b74e6961p+3, 0x1.32f1cca12a745p+3,
                            0x1.c84891698f51fp-1, -0x1.51ba54f74b9e1p+4,
                            0x1.278f8277fa44ap+4, 0x1.06260d072b40dp+4,
                            -0x1.8a63e5d6eac44p+4, 0x1.3b149e9fc56e3p+2,
                            0x1.2928780568d3bp+3, -0x1.4d2f1b585d3dfp+4,
                            0x1.3ff471d541241p+3, 0x1.fb433c2f18542p+3,
                            -0x1.c481efecb20cdp+3, -0x1.c19516511c037p-7,
                            0x1.46d2f7cbc3e2cp+2, -0x1.0f2110ced4d52p+3,
                            0x1.92004ce031f97p+1, 0x1.69f2a2d32141dp+2,
                            -0x1.fd740ba4591dcp+1, 0x1.c81649339b09cp-7,
                            0x1.13c0baf0dd59fp+0, -0x1.068ab59cb3d22p+1,
                            0x1.7aa6eeddcbffep-1, 0x1.0bfa978bc7786p+0,
                        },
                    },
                },
                { // zero 3, |z_n| = 4.2685, tau = 0.0468548, alpha = 2.1097, eta = 11.38 u
                    0x1.9e1f5267a845bp+1, 0x1.9e1f5267a845cp+1,
                    -0x1.6466d0351b28fp+1, -0x1.6466d0351b28ep+1,
                    0x1.1fc0534008766p-9,
                    {
                        { // anchor xmym
                            0x1.1584a92a59a62p-52, -0x1.98853e9fab9b9p-52,
                            0x1.df3b79d8a7f07p-52, 0x1.20dd750429b7dp+0,
                            -0x1.9227f9e6453eap+1, -0x1.d3497b09cd6d8p+1,
                            0x1.b1b3e6f3878efp+3, 0x1.4a2ec40a0832cp+0,
                            -0x1.62635a8f4915ap+4, 0x1.29b8c50b5b5e6p+4,
                            0x1.425d137f938f3p+1, -0x1.8a1facb2579b3p+5,
                            0x1.93225db3565e5p+5, 0x1.8a1e9f6a36ecdp+5,
                            -0x1.59f8c3db460f1p+6, 0x1.13ff8434876e5p+3,
                            0x1.9ada44492fc32p+5, -0x1.3dffd5250d5ddp+6,
                            0x1.f7cea5f0d5ee9p+4, 0x1.5c12313ff859cp+6,
                            -0x1.3c6ad86a1f426p+6, -0x1.6df1102bab851p+4,
                            0x1.a316b4c77e4b0p+5, -0x1.53506ff85ca05p+5,
                            0x1.27a29db104a48p+2, 0x1.97f1eacab4b95p+5,
                            -0x1.019e89d328776p+5, -0x1.0e0a098ba606dp+4,
                            0x1.4efe417d4f23cp+4, -0x1.89629362babf5p+3,
                            -0x1.64cff2e2d3944p-3, 0x1.ea7104083bc7ep+3,
                        },
                        { // anchor xmyp
                            -0x1.2c3640ddf9c8dp-52, -0x1.98853e9fab9a7p-52,
                            0x1.0798752271352p-48, 0x1.20dd750429b71p+0,
                            -0x1.9227f9e6453eep+1, -0x1.d3497b09cd6b0p+1,
                            0x1.b1b3e6f3878dcp+3, 0x1.4a2ec40a0827ap+0,
                            -0x1.62635a8f4913cp+4, 0x1.29b8c50b5b5e8p+4,
                            0x1.425d137f937cbp+1, -0x1.8a1facb2579a0p+5,
                            0x1.93225db3565e2p+5, 0x1.8a1e9f6a36ea7p+5,
                            -0x1.59f8c3db460ddp+6, 0x1.13ff84348774bp+3,
                            0x1.9ada44492fc01p+5, -0x1.3dffd5250d5d4p+6,
                            0x1.f7cea5f0d5f06p+4, 0x1.5c12313ff8583p+6,
                            -0x1.3c6ad86a1f417p+6, -0x1.6df1102bab809p+4,
                            0x1.a316b4c77e48ap+5, -0x1.53506ff85ca01p+5,
                            0x1.27a29db104ab6p+2, 0x1.97f1eacab4b7bp+5,
                            -0x1.019e89d32876bp+5, -0x1.0e0a098ba6049p+4,
                            0x1.4efe417d4f21fp+4, -0x1.89629362babf5p+3,
                            -0x1.64cff2e2d348bp-3, 0x1.ea7104083bc5ep+3,
                        },
                        { // anchor xpym
                            0x1.1584a92a59a59p-52, 0x1.526b56d14fa66p-53,
                            -0x1.56408aab303f5p-49, 0x1.20dd750429b6fp+0,
                            -0x1.9227f9e6453c2p+1, -0x1.d3497b09cd6d4p+1,
                            0x1.b1b3e6f3878d9p+3, 0x1.4a2ec40a083c1p+0,
                            -0x1.62635a8f49159p+4, 0x1.29b8c50b5b5c7p+4,
                            0x1.425d137f93a21p+1, -0x1.8a1facb2579a1p+5,
                            0x1.93225db3565bfp+5, 0x1.8a1e9f6a36ed1p+5,
                            -0x1.59f8c3db460e4p+6, 0x1.13ff843487646p+3,
                            0x1.9ada44492fc44p+5, -0x1.3dffd5250d5c5p+6,
                            0x1.f7cea5f0d5e86p+4, 0x1.5c12313ff8594p+6,
                            -0x1.3c6ad86a1f414p+6, -0x1.6df1102bab88bp+4,
                            0x1.a316b4c77e4b3p+5, -0x1.53506ff85c9dfp+5,
                            0x1.27a29db104977p+2, 0x1.97f1eacab4b88p+5,
                            -0x1.019e89d328763p+5, -0x1.0e0a098ba6083p+4,
                            0x1.4efe417d4f23cp+4, -0x1.89629362babbbp+3,
                            -0x1.64cff2e2d4125p-3, 0x1.ea7104083bc6bp+3,
                        },
                        { // anchor xpyp
                            -0x1.2c3640ddf9c78p-52, 0x1.526b56d14fa58p-53,
                            0x1.f423c17a74b2bp-51, 0x1.20dd750429b62p+0,
                            -0x1.9227f9e6453c6p+1, -0x1.d3497b09cd6acp+1,
                            0x1.b1b3e6f3878c6p+3, 0x1.4a2ec40a0830fp+0,
                            -0x1.62635a8f4913ap+4, 0x1.29b8c50b5b5c9p+4,
                            0x1.425d137f938fap+1, -0x1.8a1facb25798ep+5,
                            0x1.93225db3565bcp+5, 0x1.8a1e9f6a36eabp+5,
                            -0x1.59f8c3db460d0p+6, 0x1.13ff8434876acp+3,
                            0x1.9ada44492fc13p+5, -0x1.3dffd5250d5bcp+6,
                            0x1.f7cea5f0d5ea3p+4, 0x1.5c12313ff857cp+6,
                            -0x1.3c6ad86a1f405p+6, -0x1.6df1102bab843p+4,
                            0x1.a316b4c77e48dp+5, -0x1.53506ff85c9dbp+5,
                            0x1.27a29db1049e5p+2, 0x1.97f1eacab4b6dp+5,
                            -0x1.019e89d328759p+5, -0x1.0e0a098ba605ep+4,
                            0x1.4efe417d4f21fp+4, -0x1.89629362babbcp+3,
                            -0x1.64cff2e2d3c6dp-3, 0x1.ea7104083bc4cp+3,
                        },
                    },
                },
                { // zero 4, |z_n| = 4.9474, tau = 0.0404249, alpha = 2.1108, eta = 11.39 u
                    0x1.d94171bfa48d4p+1, 0x1.d94171bfa48d5p+1,
                    -0x1.a4c9e07279358p+1, -0x1.a4c9e07279357p+1,
                    0x1.ac63b0b7131b4p-10,
                    {
                        { // anchor xmym
                            0x1.106d57530a876p-51, -0x1.ee4ad7feeac98p-52,
                            -0x1.858d839e91adep-51, 0x1.20dd750429b89p+0,
                            -0x1.dacf1f1ea67e2p+1, -0x1.0b0182e4b359ap+2,
                            0x1.24960605224a6p+4, 0x1.66c6ac0bad6c4p+0,
                            -0x1.12094c99a1096p+5, 0x1.d8d8ffd1febdap+4,
                            0x1.1efb6bac4b6dcp+2, -0x1.653e25f129879p+6,
                            0x1.9f09014cd73aep+6, 0x1.a487a2e9a994ep+6,
                            -0x1.a345a444ac196p+7, 0x1.7d44ac98839d3p+3,
                            0x1.3c14c2f21c4b7p+7, -0x1.a32c1096afc19p+7,
                            0x1.1762395ff6176p+6, 0x1.1d01f18b87e45p+8,
                            -0x1.0ea4d7c788b81p+8, -0x1.eba6825c8ad06p+6,
                            0x1.e568831ebc686p+7, -0x1.05eb75b58fb18p+7,
                            -0x1.0596016f7d334p+5, 0x1.d453a0aa75356p+7,
                            -0x1.12569e9796f18p+7, -0x1.032be0e8028dfp+7,
                            0x1.13f58312f941bp+7, -0x1.d6cd907adc0adp+4,
                            -0x1.26aafadd8b34fp+5, 0x1.710da9eeef419p+6,
                        },
                        { // anchor xmyp
                            -0x1.0701db11f30c6p-55, -0x1.ee4ad7feeac8ep-52,
                            0x1.b49fa4e1c2462p-49, 0x1.20dd750429b7bp+0,
                            -0x1.dacf1f1ea67e6p+1, -0x1.0b0182e4b357fp+2,
                            0x1.2496060522497p+4, 0x1.66c6ac0bad5b2p+0,
                            -0x1.12094c99a107ap+5, 0x1.d8d8ffd1febdcp+4,
                            0x1.1efb6bac4b5a1p+2, -0x1.653e25f129866p+6,
                            0x1.9f09014cd73abp+6, 0x1.a487a2e9a9920p+6,
                            -0x1.a345a444ac17cp+7, 0x1.7d44ac9883b0fp+3,
                            0x1.3c14c2f21c48fp+7, -0x1.a32c1096afc0fp+7,
                            0x1.1762395ff619cp+6, 0x1.1d01f18b87e2fp+8,
                            -0x1.0ea4d7c788b76p+8, -0x1.eba6825c8acb3p+6,
                            0x1.e568831ebc65ap+7, -0x1.05eb75b58fb1ep+7,
                            -0x1.0596016f7d2cap+5, 0x1.d453a0aa7533ap+7,
                            -0x1.12569e9796f11p+7, -0x1.032be0e8028c1p+7,
                            0x1.13f58312f9405p+7, -0x1.d6cd907adc0f2p+4,
                            -0x1.26aafadd8b319p+5, 0x1.710da9eeef406p+6,
                        },
                        { // anchor xpym
                            0x1.106d57530a86bp-51, 0x1.4dc04825a29a6p-54,
                            -0x1.1e1940032573fp-48, 0x1.20dd750429b79p+0,
                            -0x1.dacf1f1ea67abp+1, -0x1.0b0182e4b3598p+2,
                            0x1.2496060522495p+4, 0x1.66c6ac0bad7b0p+0,
                            -0x1.12094c99a1094p+5, 0x1.d8d8ffd1feba2p+4,
                            0x1.1efb6bac4b814p+2, -0x1.653e25f129865p+6,
                            0x1.9f09014cd7380p+6, 0x1.a487a2e9a9951p+6,
                            -0x1.a345a444ac182p+7, 0x1.7d44ac988382fp+3,
                            0x1.3c14c2f21c4c1p+7, -0x1.a32c1096afbf1p+7,
                            0x1.1762395ff6121p+6, 0x1.1d01f18b87e3bp+8,
                            -0x1.0ea4d7c788b6cp+8, -0x1.eba6825c8ad33p+6,
                            0x1.e568831ebc680p+7, -0x1.05eb75b58faecp+7,
                            -0x1.0596016f7d3a3p+5, 0x1.d453a0aa7533cp+7,
                            -0x1.12569e9796efap+7, -0x1.032be0e8028e5p+7,
                            0x1.13f58312f9412p+7, -0x1.d6cd907adc000p+4,
                            -0x1.26aafadd8b375p+5, 0x1.710da9eeef3fep+6,
                        },
                        { // anchor xpyp
                            -0x1.0701db11f3062p-55, 0x1.4dc04825a295ap-54,
                            -0x1.317bd1e721b45p-52, 0x1.20dd750429b6ap+0,
                            -0x1.dacf1f1ea67afp+1, -0x1.0b0182e4b357dp+2,
                            0x1.2496060522486p+4, 0x1.66c6ac0bad69ep+0,
                            -0x1.12094c99a1079p+5, 0x1.d8d8ffd1feba5p+4,
                            0x1.1efb6bac4b6d8p+2, -0x1.653e25f129852p+6,
                            0x1.9f09014cd737dp+6, 0x1.a487a2e9a9923p+6,
                            -0x1.a345a444ac168p+7, 0x1.7d44ac988396cp+3,
                            0x1.3c14c2f21c499p+7, -0x1.a32c1096afbe7p+7,
                            0x1.1762395ff6147p+6, 0x1.1d01f18b87e26p+8,
                            -0x1.0ea4d7c788b61p+8, -0x1.eba6825c8ace0p+6,
                            0x1.e568831ebc654p+7, -0x1.05eb75b58faf2p+7,
                            -0x1.0596016f7d33ap+5, 0x1.d453a0aa75320p+7,
                            -0x1.12569e9796ef3p+7, -0x1.032be0e8028c7p+7,
                            0x1.13f58312f93fdp+7, -0x1.d6cd907adc045p+4,
                            -0x1.26aafadd8b33fp+5, 0x1.710da9eeef3ebp+6,
                        },
                    },
                },
                { // zero 5, |z_n| = 5.5446, tau = 0.036071, alpha = 2.1114, eta = 11.39 u
                    0x1.06ca76356512bp+2, 0x1.06ca76356512cp+2,
                    -0x1.dcebe33c4bdf4p+1, -0x1.dcebe33c4bdf3p+1,
                    0x1.55147f45c44bep-10,
                    {
                        { // anchor xmym
                            0x1.b25764272b243p-52, -0x1.f11ebfa493224p-53,
                            -0x1.acaa459e558a1p-50, 0x1.20dd750429b82p+0,
                            -0x1.0d12f8a348556p+2, -0x1.288719e2ebd15p+2,
                            0x1.70484a3d8e2fap+4, 0x1.7cce3f2ab7f00p+0,
                            -0x1.7f67828fe2143p+5, 0x1.5126cf5147425p+5,
                            0x1.aca2ab1b48382p+2, -0x1.1a7b89cf9eb48p+7,
                            0x1.6c748f6a2fd38p+7, 0x1.772cbf9513370p+7,
                            -0x1.9f6549299fff7p+8, 0x1.c8c2be759df1cp+3,
                            0x1.6f8fe0c50321bp+8, -0x1.c07c27acbaa58p+8,
                            0x1.008618be63138p+7, 0x1.63325f7826348p+9,
                            -0x1.621c6468fc5a5p+9, -0x1.8e1d2ca6ce514p+8,
                            0x1.838e8bf2aadb1p+9, -0x1.37b9dbd587d2ep+8,
                            -0x1.b5aae8c77e76dp+7, 0x1.7c82bb9b94ec7p+9,
                            -0x1.a13e6d866a9f5p+8, -0x1.171d00f1ad8e8p+9,
                            0x1.1e9262ed43719p+9, -0x1.adf06d6d96fa5p+1,
                            -0x1.007bea6febdd5p+8, 0x1.6900c3e9e7f85p+8,
                        },
                        { // anchor xmyp
                            -0x1.1ec70bc25095fp-53, -0x1.f11ebfa49321ap-53,
                            0x1.7ab910f6acdb7p-49, 0x1.20dd750429b71p+0,
                            -0x1.0d12f8a348558p+2, -0x1.288719e2ebcf3p+2,
                            0x1.70484a3d8e2e4p+4, 0x1.7cce3f2ab7d81p+0,
                            -0x1.7f67828fe2117p+5, 0x1.5126cf5147427p+5,
                            0x1.aca2ab1b4814fp+2, -0x1.1a7b89cf9eb37p+7,
                            0x1.6c748f6a2fd37p+7, 0x1.772cbf9513343p+7,
                            -0x1.9f6549299ffdbp+8, 0x1.c8c2be759e1fbp+3,
                            0x1.6f8fe0c5031e9p+8, -0x1.c07c27acbaa4fp+8,
                            0x1.008618be63176p+7, 0x1.63325f782632cp+9,
                            -0x1.621c6468fc598p+9, -0x1.8e1d2ca6ce4d1p+8,
                            0x1.838e8bf2aad8dp+9, -0x1.37b9dbd587d43p+8,
                            -0x1.b5aae8c77e6fcp+7, 0x1.7c82bb9b94eb1p+9,
                            -0x1.a13e6d866a9f5p+8, -0x1.171d00f1ad8c8p+9,
                            0x1.1e9262ed43704p+9, -0x1.adf06d6d97eacp+1,
                            -0x1.007bea6febdafp+8, 0x1.6900c3e9e7f77p+8,
                        },
                        { // anchor xpym
                            0x1.b25764272b1e5p-52, 0x1.c5733a1f2ea55p-51,
                            -0x1.42a8415713047p-47, 0x1.20dd750429b5dp+0,
                            -0x1.0d12f8a348511p+2, -0x1.288719e2ebd11p+2,
                            0x1.70484a3d8e2cap+4, 0x1.7cce3f2ab81a3p+0,
                            -0x1.7f67828fe213fp+5, 0x1.5126cf51473cdp+5,
                            0x1.aca2ab1b487c7p+2, -0x1.1a7b89cf9eb25p+7,
                            0x1.6c748f6a2fcdep+7, 0x1.772cbf9513373p+7,
                            -0x1.9f6549299ffc9p+8, 0x1.c8c2be759d81ap+3,
                            0x1.6f8fe0c50322dp+8, -0x1.c07c27acba9f4p+8,
                            0x1.008618be6305bp+7, 0x1.63325f7826329p+9,
                            -0x1.621c6468fc563p+9, -0x1.8e1d2ca6ce549p+8,
                            0x1.838e8bf2aad9dp+9, -0x1.37b9dbd587ca0p+8,
                            -0x1.b5aae8c77e816p+7, 0x1.7c82bb9b94e8ep+9,
                            -0x1.a13e6d866a978p+8, -0x1.171d00f1ad8e8p+9,
                            0x1.1e9262ed436fbp+9, -0x1.adf06d6d94557p+1,
                            -0x1.007bea6febdf0p+8, 0x1.6900c3e9e7f39p+8,
                        },
                        { // anchor xpyp
                            -0x1.1ec70bc250987p-53, 0x1.c5733a1f2ea36p-51,
                            -0x1.5cc968cb3a38fp-48, 0x1.20dd750429b4cp+0,
                            -0x1.0d12f8a348513p+2, -0x1.288719e2ebceep+2,
                            0x1.70484a3d8e2b5p+4, 0x1.7cce3f2ab8023p+0,
                            -0x1.7f67828fe2113p+5, 0x1.5126cf51473cfp+5,
                            0x1.aca2ab1b48594p+2, -0x1.1a7b89cf9eb14p+7,
                            0x1.6c748f6a2fcdcp+7, 0x1.772cbf9513346p+7,
                            -0x1.9f6549299ffadp+8, 0x1.c8c2be759daf9p+3,
                            0x1.6f8fe0c5031fbp+8, -0x1.c07c27acba9ebp+8,
                            0x1.008618be63099p+7, 0x1.63325f782630dp+9,
                            -0x1.621c6468fc555p+9, -0x1.8e1d2ca6ce507p+8,
                            0x1.838e8bf2aad79p+9, -0x1.37b9dbd587cb4p+8,
                            -0x1.b5aae8c77e7a5p+7, 0x1.7c82bb9b94e79p+9,
                            -0x1.a13e6d866a978p+8, -0x1.171d00f1ad8c9p+9,
                            0x1.1e9262ed436e6p+9, -0x1.adf06d6d9545ep+1,
                            -0x1.007bea6febdcap+8, 0x1.6900c3e9e7f2bp+8,
                        },
                    },
                },
                { // zero 6, |z_n| = 6.0839, tau = 0.0328739, alpha = 2.1117, eta = 11.39 u
                    0x1.1e8425f7e0f8cp+2, 0x1.1e8425f7e0f8dp+2,
                    -0x1.07a81a8382ae8p+2, -0x1.07a81a8382ae7p+2,
                    0x1.1b4c1b3bd46b5p-10,
                    {
                        { // anchor xmym
                            0x1.16b17e326a090p-51, -0x1.07da2d0bcc750p-50,
                            0x1.cf254b73bf134p-49, 0x1.20dd750429ba4p+0,
                            -0x1.29813276f4457p+2, -0x1.434c85c2b2329p+2,
                            0x1.bbf536d38a29ep+4, 0x1.8eba20eb2ae0ep+0,
                            -0x1.f7f37d15a9cc4p+5, 0x1.c18d1e5846b61p+5,
                            0x1.23688bc36d8cdp+3, -0x1.9a1c7ed57fde4p+7,
                            0x1.20fe4fcf78cbbp+8, 0x1.2bc5ccfa67ccfp+8,
                            -0x1.6a8c1863f1893p+9, 0x1.e9d6bee2dbcbep+3,
                            0x1.69c1ea3088772p+9, -0x1.a36da33cfdc41p+9,
                            0x1.a297d17933ca8p+7, 0x1.755f85a26e92fp+10,
                            -0x1.86a83592c8a38p+10, -0x1.f2803a983e94dp+9,
                            0x1.ef29ba28470b5p+10, -0x1.3b3b86998f8c8p+9,
                            -0x1.8842200cfb85fp+9, 0x1.f32070c40e001p+10,
                            -0x1.0170348311fcdp+10, -0x1.bbd3348b3c653p+10,
                            0x1.c5dc88260d7c8p+10, 0x1.e8477240a6e89p+7,
                            -0x1.0e1d42a7fbd21p+10, 0x1.100b26dde098fp+10,
                        },
                        { // anchor xmyp
                            -0x1.2b096bd5e9693p-51, -0x1.07da2d0bcc72dp-50,
                            0x1.b715d89fa1f4cp-47, 0x1.20dd750429b7fp+0,
                            -0x1.29813276f445cp+2, -0x1.434c85c2b22d6p+2,
                            0x1.bbf536d38a266p+4, 0x1.8eba20eb2aa1ep+0,
                            -0x1.f7f37d15a9c44p+5, 0x1.c18d1e5846b67p+5,
                            0x1.23688bc36d54ap+3, -0x1.9a1c7ed57fdadp+7,
                            0x1.20fe4fcf78cbap+8, 0x1.2bc5ccfa67c80p+8,
                            -0x1.6a8c1863f185ep+9, 0x1.e9d6bee2dc80cp+3,
                            0x1.69c1ea3088709p+9, -0x1.a36da33cfdc32p+9,
                            0x1.a297d17933ddfp+7, 0x1.755f85a26e8f2p+10,
                            -0x1.86a83592c8a1dp+10, -0x1.f2803a983e8a2p+9,
                            0x1.ef29ba2847058p+10, -0x1.3b3b86998f911p+9,
                            -0x1.8842200cfb7abp+9, 0x1.f32070c40dfcdp+10,
                            -0x1.0170348311fdbp+10, -0x1.bbd3348b3c5f0p+10,
                            0x1.c5dc88260d788p+10, 0x1.e8477240a6c8ep+7,
                            -0x1.0e1d42a7fbcd6p+10, 0x1.100b26dde0983p+10,
                        },
                        { // anchor xpym
                            0x1.16b17e326a088p-51, 0x1.90347f85d43fep-54,
                            -0x1.6b6fbf3408fc1p-48, 0x1.20dd750429b7cp+0,
                            -0x1.29813276f4404p+2, -0x1.434c85c2b2324p+2,
                            0x1.bbf536d38a25fp+4, 0x1.8eba20eb2b191p+0,
                            -0x1.f7f37d15a9cbep+5, 0x1.c18d1e5846ae1p+5,
                            0x1.23688bc36dc30p+3, -0x1.9a1c7ed57fdabp+7,
                            0x1.20fe4fcf78c6cp+8, 0x1.2bc5ccfa67cd1p+8,
                            -0x1.6a8c1863f1866p+9, 0x1.e9d6bee2dafa3p+3,
                            0x1.69c1ea3088781p+9, -0x1.a36da33cfdbd8p+9,
                            0x1.a297d17933abfp+7, 0x1.755f85a26e908p+10,
                            -0x1.86a83592c89e3p+10, -0x1.f2803a983e983p+9,
                            0x1.ef29ba2847091p+10, -0x1.3b3b86998f80cp+9,
                            -0x1.8842200cfb8c8p+9, 0x1.f32070c40dfa7p+10,
                            -0x1.0170348311f6ap+10, -0x1.bbd3348b3c646p+10,
                            0x1.c5dc88260d788p+10, 0x1.e8477240a7087p+7,
                            -0x1.0e1d42a7fbd2cp+10, 0x1.100b26dde0944p+10,
                        },
                        { // anchor xpyp
                            -0x1.2b096bd5e964ap-51, 0x1.90347f85d43d6p-54,
                            0x1.1b294c515b635p-48, 0x1.20dd750429b56p+0,
                            -0x1.29813276f4408p+2, -0x1.434c85c2b22d1p+2,
                            0x1.bbf536d38a227p+4, 0x1.8eba20eb2ada1p+0,
                            -0x1.f7f37d15a9c3ep+5, 0x1.c18d1e5846ae7p+5,
                            0x1.23688bc36d8adp+3, -0x1.9a1c7ed57fd75p+7,
                            0x1.20fe4fcf78c6ap+8, 0x1.2bc5ccfa67c81p+8,
                            -0x1.6a8c1863f1831p+9, 0x1.e9d6bee2dbaf1p+3,
                            0x1.69c1ea3088718p+9, -0x1.a36da33cfdbc9p+9,
                            0x1.a297d17933bf7p+7, 0x1.755f85a26e8cbp+10,
                            -0x1.86a83592c89c7p+10, -0x1.f2803a983e8d9p+9,
                            0x1.ef29ba2847033p+10, -0x1.3b3b86998f856p+9,
                            -0x1.8842200cfb814p+9, 0x1.f32070c40df73p+10,
                            -0x1.0170348311f77p+10, -0x1.bbd3348b3c5e3p+10,
                            0x1.c5dc88260d749p+10, 0x1.e8477240a6e8cp+7,
                            -0x1.0e1d42a7fbce1p+10, 0x1.100b26dde0939p+10,
                        },
                    },
                },
                { // zero 7, |z_n| = 6.5793, tau = 0.0303985, alpha = 2.1119, eta = 11.39 u
                    0x1.34621cb7730e7p+2, 0x1.34621cb7730e8p+2,
                    -0x1.1eb5949f736fdp+2, -0x1.1eb5949f736fcp+2,
                    0x1.e47aab66a0bf2p-11,
                    {
                        { // anchor xmym
                            0x1.7e1e966bce431p-51, -0x1.1af76d6b99a7ep-51,
                            -0x1.1ecc3c2d306e4p-49, 0x1.20dd750429b9dp+0,
                            -0x1.43844d0256023p+2, -0x1.5bf925f706bb0p+2,
                            0x1.03cf710f0fd98p+5, 0x1.9dd5596064251p+0,
                            -0x1.3d58d33bca689p+6, 0x1.1e45b7b23b575p+6,
                            0x1.75ae584589b8fp+3, -0x1.18c0b3acd2f56p+8,
                            0x1.aaee6c86fc583p+8, 0x1.bc8430c343f40p+8,
                            -0x1.2203a9441f36ep+10, 0x1.d647a3210b929p+3,
                            0x1.3e8f4398b4299p+10, -0x1.650313f49f467p+10,
                            0x1.3b0404156e02cp+8, 0x1.5d4c2cf6c03cep+11,
                            -0x1.7ec32b7111e4fp+11, -0x1.09a2b6a0c7c2cp+11,
                            0x1.1044f0073c57bp+12, -0x1.1d1ba610aab59p+10,
                            -0x1.0b142683f0f73p+11, 0x1.1aaa12970b604p+12,
                            -0x1.136a602596a25p+11, -0x1.229c6dd33f495p+12,
                            0x1.2bda6320a158ap+12, 0x1.1e22350c1896ep+10,
                            -0x1.b20609ff7d536p+11, 0x1.57cabd43b11b4p+11,
                        },
                        { // anchor xmyp
                            -0x1.8738a7390a5c2p-52, -0x1.1af76d6b99a65p-51,
                            0x1.144616ebba9c6p-47, 0x1.20dd750429b75p+0,
                            -0x1.43844d0256028p+2, -0x1.5bf925f706b4fp+2,
                            0x1.03cf710f0fd75p+5, 0x1.9dd5596063d5cp+0,
                            -0x1.3d58d33bca631p+6, 0x1.1e45b7b23b579p+6,
                            0x1.75ae584589659p+3, -0x1.18c0b3acd2f2ep+8,
                            0x1.aaee6c86fc581p+8, 0x1.bc8430c343ec1p+8,
                            -0x1.2203a9441f341p+10, 0x1.d647a3210cd12p+3,
                            0x1.3e8f4398b4237p+10, -0x1.650313f49f45bp+10,
                            0x1.3b0404156e178p+8, 0x1.5d4c2cf6c0393p+11,
                            -0x1.7ec32b7111e37p+11, -0x1.09a2b6a0c7bcfp+11,
                            0x1.1044f0073c546p+12, -0x1.1d1ba610aabbdp+10,
                            -0x1.0b142683f0efdp+11, 0x1.1aaa12970b5e8p+12,
                            -0x1.136a602596a44p+11, -0x1.229c6dd33f454p+12,
                            0x1.2bda6320a1562p+12, 0x1.1e22350c188a2p+10,
                            -0x1.b20609ff7d4c1p+11, 0x1.57cabd43b11b2p+11,
                        },
                        { // anchor xpym
                            0x1.7e1e966bce3f6p-51, 0x1.26c37c9cb9c91p-51,
                            -0x1.8b375c0da21acp-47, 0x1.20dd750429b72p+0,
                            -0x1.43844d0255fc2p+2, -0x1.5bf925f706babp+2,
                            0x1.03cf710f0fd71p+5, 0x1.9dd55960646cap+0,
                            -0x1.3d58d33bca686p+6, 0x1.1e45b7b23b51ep+6,
                            0x1.75ae58458a090p+3, -0x1.18c0b3acd2f2cp+8,
                            0x1.aaee6c86fc504p+8, 0x1.bc8430c343f42p+8,
                            -0x1.2203a9441f346p+10, 0x1.d647a3210a2d9p+3,
                            0x1.3e8f4398b42a4p+10, -0x1.650313f49f404p+10,
                            0x1.3b0404156de4ep+8, 0x1.5d4c2cf6c03a5p+11,
                            -0x1.7ec32b7111df1p+11, -0x1.09a2b6a0c7c45p+11,
                            0x1.1044f0073c562p+12, -0x1.1d1ba610aaa85p+10,
                            -0x1.0b142683f0fabp+11, 0x1.1aaa12970b5c9p+12,
                            -0x1.136a6025969a1p+11, -0x1.229c6dd33f486p+12,
                            0x1.2bda6320a1558p+12, 0x1.1e22350c18a0fp+10,
                            -0x1.b20609ff7d537p+11, 0x1.57cabd43b113fp+11,
                        },
                        { // anchor xpyp
                            -0x1.8738a7390a589p-52, 0x1.26c37c9cb9c59p-51,
                            -0x1.79f1b0b4db188p-50, 0x1.20dd750429b49p+0,
                            -0x1.43844d0255fc7p+2, -0x1.5bf925f706b4ap+2,
                            0x1.03cf710f0fd4dp+5, 0x1.9dd55960641d5p+0,
                            -0x1.3d58d33bca62ep+6, 0x1.1e45b7b23b521p+6,
                            0x1.75ae584589b5ap+3, -0x1.18c0b3acd2f04p+8,
                            0x1.aaee6c86fc502p+8, 0x1.bc8430c343ec3p+8,
                            -0x1.2203a9441f319p+10, 0x1.d647a3210b6c2p+3,
                            0x1.3e8f4398b4242p+10, -0x1.650313f49f3f9p+10,
                            0x1.3b0404156df9ap+8, 0x1.5d4c2cf6c0369p+11,
                            -0x1.7ec32b7111dd9p+11, -0x1.09a2b6a0c7be7p+11,
                            0x1.1044f0073c52dp+12, -0x1.1d1ba610aaae9p+10,
                            -0x1.0b142683f0f35p+11, 0x1.1aaa12970b5adp+12,
                            -0x1.136a6025969c1p+11, -0x1.229c6dd33f444p+12,
                            0x1.2bda6320a152fp+12, 0x1.1e22350c18943p+10,
                            -0x1.b20609ff7d4c2p+11, 0x1.57cabd43b113dp+11,
                        },
                    },
                },
                { // zero 8, |z_n| = 7.0400, tau = 0.0284089, alpha = 2.1120, eta = 11.39 u
                    0x1.48c5b5cddc381p+2, 0x1.48c5b5cddc382p+2,
                    -0x1.3415689f296a6p+2, -0x1.3415689f296a5p+2,
                    0x1.a722aaaafb341p-11,
                    {
                        { // anchor xmym
                            0x1.39ee08368ce2bp-51, -0x1.3e3291c26bcbdp-52,
                            -0x1.a76767b9a92cdp-49, 0x1.20dd750429b91p+0,
                            -0x1.5ba299012b59fp+2, -0x1.72fad57d3f3a7p+2,
                            0x1.29a334b1d540cp+5, 0x1.aae3c24036302p+0,
                            -0x1.836fd279e9b4dp+6, 0x1.60a72a0e0288bp+6,
                            0x1.cc5fbe1ab4375p+3, -0x1.7054f1c908ce3p+8,
                            0x1.2b56b28786af7p+9, 0x1.3834c515d644ep+9,
                            -0x1.b366fd9ca45e8p+10, 0x1.857e7653df6a6p+3,
                            0x1.030bf452842aep+11, -0x1.1b75ac69cc8d8p+11,
                            0x1.bf13bbee979b8p+8, 0x1.2c2e61c0e62acp+12,
                            -0x1.579da62133157p+12, -0x1.fa532a5d63518p+11,
                            0x1.0cb817dc69c6cp+13, -0x1.db77df4287969p+10,
                            -0x1.3381991e3529cp+12, 0x1.1f93811911c4dp+13,
                            -0x1.099acb6733b11p+12, -0x1.4c00192fe0065p+13,
                            0x1.5bbc76bf9775cp+13, 0x1.bcf0fc44b87d1p+11,
                            -0x1.23de1d52d5525p+13, 0x1.7e8955f28c70ap+12,
                        },
                        { // anchor xmyp
                            -0x1.07cce1d1c68ccp-51, -0x1.3e3291c26bc95p-52,
                            0x1.0920fb8ed4ebcp-47, 0x1.20dd750429b66p+0,
                            -0x1.5ba299012b5a4p+2, -0x1.72fad57d3f338p+2,
                            0x1.29a334b1d53e0p+5, 0x1.aae3c24035cf5p+0,
                            -0x1.836fd279e9ad9p+6, 0x1.60a72a0e0288fp+6,
                            0x1.cc5fbe1ab3c24p+3, -0x1.7054f1c908cabp+8,
                            0x1.2b56b28786af6p+9, 0x1.3834c515d63efp+9,
                            -0x1.b366fd9ca45a2p+10, 0x1.857e7653e1708p+3,
                            0x1.030bf4528425ap+11, -0x1.1b75ac69cc8d0p+11,
                            0x1.bf13bbee97c31p+8, 0x1.2c2e61c0e6276p+12,
                            -0x1.579da62133142p+12, -0x1.fa532a5d6345fp+11,
                            0x1.0cb817dc69c37p+13, -0x1.db77df4287a4fp+10,
                            -0x1.3381991e35215p+12, 0x1.1f93811911c32p+13,
                            -0x1.099acb6733b42p+12, -0x1.4c00192fe0019p+13,
                            0x1.5bbc76bf9772fp+13, 0x1.bcf0fc44b86bfp+11,
                            -0x1.23de1d52d54d7p+13, 0x1.7e8955f28c719p+12,
                        },
                        { // anchor xpym
                            0x1.39ee08368cde5p-51, 0x1.a2a1a1271d896p-51,
                            -0x1.c57c72ef95a1bp-47, 0x1.20dd750429b63p+0,
                            -0x1.5ba299012b530p+2, -0x1.72fad57d3f3a2p+2,
                            0x1.29a334b1d53dbp+5, 0x1.aae3c24036885p+0,
                            -0x1.836fd279e9b48p+6, 0x1.60a72a0e02817p+6,
                            0x1.cc5fbe1ab4a79p+3, -0x1.7054f1c908ca9p+8,
                            0x1.2b56b28786a98p+9, 0x1.3834c515d644fp+9,
                            -0x1.b366fd9ca45a8p+10, 0x1.857e7653dd338p+3,
                            0x1.030bf452842b6p+11, -0x1.1b75ac69cc883p+11,
                            0x1.bf13bbee9765dp+8, 0x1.2c2e61c0e6284p+12,
                            -0x1.579da621330fap+12, -0x1.fa532a5d63541p+11,
                            0x1.0cb817dc69c50p+13, -0x1.db77df42877b9p+10,
                            -0x1.3381991e352d2p+12, 0x1.1f93811911c09p+13,
                            -0x1.099acb6733a79p+12, -0x1.4c00192fe004cp+13,
                            0x1.5bbc76bf97718p+13, 0x1.bcf0fc44b8884p+11,
                            -0x1.23de1d52d551ep+13, 0x1.7e8955f28c66ep+12,
                        },
                        { // anchor xpyp
                            -0x1.07cce1d1c68b5p-51, 0x1.a2a1a1271d853p-51,
                            -0x1.4a0675c959ac1p-49, 0x1.20dd750429b38p+0,
                            -0x1.5ba299012b535p+2, -0x1.72fad57d3f333p+2,
                            0x1.29a334b1d53afp+5, 0x1.aae3c24036277p+0,
                            -0x1.836fd279e9ad5p+6, 0x1.60a72a0e0281cp+6,
                            0x1.cc5fbe1ab4328p+3, -0x1.7054f1c908c71p+8,
                            0x1.2b56b28786a97p+9, 0x1.3834c515d63f0p+9,
                            -0x1.b366fd9ca4561p+10, 0x1.857e7653df399p+3,
                            0x1.030bf45284261p+11, -0x1.1b75ac69cc87bp+11,
                            0x1.bf13bbee978d6p+8, 0x1.2c2e61c0e624ep+12,
                            -0x1.579da621330e6p+12, -0x1.fa532a5d63488p+11,
                            0x1.0cb817dc69c1ap+13, -0x1.db77df42878a0p+10,
                            -0x1.3381991e3524bp+12, 0x1.1f93811911beep+13,
                            -0x1.099acb6733aaap+12, -0x1.4c00192fe0000p+13,
                            0x1.5bbc76bf976ebp+13, 0x1.bcf0fc44b8773p+11,
                            -0x1.23de1d52d54d0p+13, 0x1.7e8955f28c67dp+12,
                        },
                    },
                },
                { // zero 9, |z_n| = 7.4726, tau = 0.0267646, alpha = 2.1121, eta = 11.39 u
                    0x1.5bf2f4f5a8d21p+2, 0x1.5bf2f4f5a8d22p+2,
                    -0x1.481917c474953p+2, -0x1.481917c474952p+2,
                    0x1.779218c8237c5p-11,
                    {
                        { // anchor xmym
                            0x1.3a8f5d1b5c962p-53, -0x1.299e47a88840cp-52,
                            0x1.4f5514147ca19p-50, 0x1.20dd750429b80p+0,
                            -0x1.72380e5f5579fp+2, -0x1.889e535d8872cp+2,
                            0x1.4f7639d1004fcp+5, 0x1.b662dda32ae77p+0,
                            -0x1.cdf0abcf30d54p+6, 0x1.a7968de483cb1p+6,
                            0x1.1375f576d73f1p+4, -0x1.d3cad757f08b4p+8,
                            0x1.9356ed3949874p+9, 0x1.a4ef6b08d1966p+9,
                            -0x1.3762717b25f0cp+11, 0x1.e024eb37934f8p+2,
                            0x1.8cd0d6d379fc3p+11, -0x1.aaa9ea78cd6b4p+11,
                            0x1.2f653433ddc5fp+9, 0x1.e3aa87cf5f1c5p+12,
                            -0x1.2069820dd0921p+13, -0x1.bc5c017014a54p+12,
                            0x1.e8bcfeb4dd296p+13, -0x1.74b3105fd1afbp+11,
                            -0x1.3b2cc3db7f25ep+13, 0x1.0d85e335864cep+14,
                            -0x1.d96f3bf5c766ap+12, -0x1.568e07d334dacp+14,
                            0x1.6d4ff473d9987p+14, 0x1.19bc79ee241d3p+13,
                            -0x1.5951f5abc94fdp+14, 0x1.828748edfacfep+13,
                        },
                        { // anchor xmyp
                            -0x1.f31712c17c479p-51, -0x1.299e47a888395p-52,
                            0x1.b288f5e018030p-47, 0x1.20dd750429b52p+0,
                            -0x1.72380e5f557a4p+2, -0x1.889e535d886aep+2,
                            0x1.4f7639d1004c8p+5, 0x1.b662dda32a73fp+0,
                            -0x1.cdf0abcf30cc2p+6, 0x1.a7968de483cb6p+6,
                            0x1.1375f576d6f02p+4, -0x1.d3cad757f0868p+8,
                            0x1.9356ed3949873p+9, 0x1.a4ef6b08d18ddp+9,
                            -0x1.3762717b25ed6p+11, 0x1.e024eb379982cp+2,
                            0x1.8cd0d6d379f3bp+11, -0x1.aaa9ea78cd6a9p+11,
                            0x1.2f653433dde8bp+9, 0x1.e3aa87cf5f16bp+12,
                            -0x1.2069820dd0911p+13, -0x1.bc5c0170149acp+12,
                            0x1.e8bcfeb4dd231p+13, -0x1.74b3105fd1be7p+11,
                            -0x1.3b2cc3db7f1d3p+13, 0x1.0d85e335864b6p+14,
                            -0x1.d96f3bf5c76e5p+12, -0x1.568e07d334d5cp+14,
                            0x1.6d4ff473d995ap+14, 0x1.19bc79ee24131p+13,
                            -0x1.5951f5abc94a0p+14, 0x1.828748edfad1ep+13,
                        },
                        { // anchor xpym
                            0x1.3a8f5d1b5c8d2p-53, 0x1.acebc6340f4c9p-51,
                            -0x1.484d6bdcc5e1dp-47, 0x1.20dd750429b4fp+0,
                            -0x1.72380e5f55721p+2, -0x1.889e535d88726p+2,
                            0x1.4f7639d1004c3p+5, 0x1.b662dda32b515p+0,
                            -0x1.cdf0abcf30d4ep+6, 0x1.a7968de483c1fp+6,
                            0x1.1375f576d78abp+4, -0x1.d3cad757f0865p+8,
                            0x1.9356ed39497ecp+9, 0x1.a4ef6b08d1966p+9,
                            -0x1.3762717b25edap+11, 0x1.e024eb378ca4ep+2,
                            0x1.8cd0d6d379fcep+11, -0x1.aaa9ea78cd62cp+11,
                            0x1.2f653433dd98ep+9, 0x1.e3aa87cf5f17fp+12,
                            -0x1.2069820dd08cdp+13, -0x1.bc5c017014a74p+12,
                            0x1.e8bcfeb4dd25bp+13, -0x1.74b3105fd1966p+11,
                            -0x1.3b2cc3db7f28ep+13, 0x1.0d85e33586489p+14,
                            -0x1.d96f3bf5c752ap+12, -0x1.568e07d334d8dp+14,
                            0x1.6d4ff473d9936p+14, 0x1.19bc79ee2422dp+13,
                            -0x1.5951f5abc94edp+14, 0x1.828748edfac44p+13,
                        },
                        { // anchor xpyp
                            -0x1.f31712c17c43bp-51, 0x1.acebc6340f4a8p-51,
                            0x1.01439e030a32ap-49, 0x1.20dd750429b21p+0,
                            -0x1.72380e5f55726p+2, -0x1.889e535d886a9p+2,
                            0x1.4f7639d10048ep+5, 0x1.b662dda32addep+0,
                            -0x1.cdf0abcf30cbcp+6, 0x1.a7968de483c24p+6,
                            0x1.1375f576d73bcp+4, -0x1.d3cad757f0819p+8,
                            0x1.9356ed39497ebp+9, 0x1.a4ef6b08d18dep+9,
                            -0x1.3762717b25ea5p+11, 0x1.e024eb3792d82p+2,
                            0x1.8cd0d6d379f46p+11, -0x1.aaa9ea78cd621p+11,
                            0x1.2f653433ddbbap+9, 0x1.e3aa87cf5f125p+12,
                            -0x1.2069820dd08bdp+13, -0x1.bc5c0170149ccp+12,
                            0x1.e8bcfeb4dd1f6p+13, -0x1.74b3105fd1a53p+11,
                            -0x1.3b2cc3db7f203p+13, 0x1.0d85e33586470p+14,
                            -0x1.d96f3bf5c75a6p+12, -0x1.568e07d334d3dp+14,
                            0x1.6d4ff473d9909p+14, 0x1.19bc79ee2418bp+13,
                            -0x1.5951f5abc9490p+14, 0x1.828748edfac65p+13,
                        },
                    },
                },
                { // zero 10, |z_n| = 7.8814, tau = 0.0253761, alpha = 2.1122, eta = 11.39 u
                    0x1.6e1b9ac700648p+2, 0x1.6e1b9ac700649p+2,
                    -0x1.5afb4ea79016cp+2, -0x1.5afb4ea79016bp+2,
                    0x1.519cfc706ae32p-11,
                    {
                        { // anchor xmym
                            0x1.117b13c7b3bbbp-50, -0x1.0c10fa4c4da56p-52,
                            -0x1.2c45ed47ef05cp-47, 0x1.20dd750429ba7p+0,
                            -0x1.8786e209e82c8p+2, -0x1.9d1bc34d10180p+2,
                            0x1.7548b4bcff9e0p+5, 0x1.c0a7f08a45338p+0,
                            -0x1.0e4d6b7ee1b30p+7, 0x1.f2ce8f849c28ep+6,
                            0x1.4271814bef44cp+4, -0x1.2191256536c68p+9,
                            0x1.074ed372a3975p+10, 0x1.12c9ba9886c95p+10,
                            -0x1.acb91ab403364p+11, 0x1.eff61cb92f9f8p-2,
                            0x1.22100156b5d66p+12, -0x1.33c3366da7761p+12,
                            0x1.8d94e9dad98d1p+9, 0x1.725851d6ce459p+13,
                            -0x1.cb01d465fc28dp+13, -0x1.6dfc90611a868p+13,
                            0x1.a0d3c9cb46632p+14, -0x1.1649974d86de3p+12,
                            -0x1.284a3b2d44cdbp+14, 0x1.d9793464201c9p+14,
                            -0x1.8c9583b08234ep+13, -0x1.468f0def3f1b7p+15,
                            0x1.631c3028c9987p+15, 0x1.388320e760576p+14,
                            -0x1.7296007052b4fp+15, 0x1.6a1fde4d7bbb3p+14,
                        },
                        { // anchor xmyp
                            -0x1.ec4c278ebfa7cp-55, -0x1.0c10fa4c4da85p-52,
                            0x1.c357581484378p-49, 0x1.20dd750429b76p+0,
                            -0x1.8786e209e82cdp+2, -0x1.9d1bc34d100f4p+2,
                            0x1.7548b4bcff9a2p+5, 0x1.c0a7f08a44ac6p+0,
                            -0x1.0e4d6b7ee1ad5p+7, 0x1.f2ce8f849c294p+6,
                            0x1.4271814beeddbp+4, -0x1.2191256536c36p+9,
                            0x1.074ed372a3975p+10, 0x1.12c9ba9886c38p+10,
                            -0x1.acb91ab403317p+11, 0x1.eff61cb9c0a78p-2,
                            0x1.22100156b5cfep+12, -0x1.33c3366da775ap+12,
                            0x1.8d94e9dad9c64p+9, 0x1.725851d6ce411p+13,
                            -0x1.cb01d465fc275p+13, -0x1.6dfc90611a7d9p+13,
                            0x1.a0d3c9cb465d9p+14, -0x1.1649974d86ec1p+12,
                            -0x1.284a3b2d44c56p+14, 0x1.d9793464201a1p+14,
                            -0x1.8c9583b0823d6p+13, -0x1.468f0def3f169p+15,
                            0x1.631c3028c995cp+15, 0x1.388320e7604c8p+14,
                            -0x1.7296007052aebp+15, 0x1.6a1fde4d7bbe4p+14,
                        },
                        { // anchor xpym
                            0x1.117b13c7b3b7dp-50, 0x1.bbb26ce22c9f0p-51,
                            -0x1.59e667a8eb96fp-46, 0x1.20dd750429b74p+0,
                            -0x1.8786e209e823cp+2, -0x1.9d1bc34d1017bp+2,
                            0x1.7548b4bcff99dp+5, 0x1.c0a7f08a45b03p+0,
                            -0x1.0e4d6b7ee1b2dp+7, 0x1.f2ce8f849c1d9p+6,
                            0x1.4271814befa78p+4, -0x1.2191256536c34p+9,
                            0x1.074ed372a3917p+10, 0x1.12c9ba9886c95p+10,
                            -0x1.acb91ab40331bp+11, 0x1.eff61cb895bdep-2,
                            0x1.22100156b5d6dp+12, -0x1.33c3366da76f9p+12,
                            0x1.8d94e9dad9456p+9, 0x1.725851d6ce420p+13,
                            -0x1.cb01d465fc1fdp+13, -0x1.6dfc90611a880p+13,
                            0x1.a0d3c9cb465fap+14, -0x1.1649974d86c80p+12,
                            -0x1.284a3b2d44d03p+14, 0x1.d979346420144p+14,
                            -0x1.8c9583b082217p+13, -0x1.468f0def3f195p+15,
                            0x1.631c3028c9930p+15, 0x1.388320e7605cbp+14,
                            -0x1.7296007052b37p+15, 0x1.6a1fde4d7baebp+14,
                        },
                        { // anchor xpyp
                            -0x1.ec4c278ebfbc8p-55, 0x1.bbb26ce22c977p-51,
                            -0x1.16b10c04c71a9p-47, 0x1.20dd750429b43p+0,
                            -0x1.8786e209e8241p+2, -0x1.9d1bc34d100efp+2,
                            0x1.7548b4bcff95ep+5, 0x1.c0a7f08a45291p+0,
                            -0x1.0e4d6b7ee1ad2p+7, 0x1.f2ce8f849c1dfp+6,
                            0x1.4271814bef407p+4, -0x1.2191256536c03p+9,
                            0x1.074ed372a3917p+10, 0x1.12c9ba9886c38p+10,
                            -0x1.acb91ab4032cep+11, 0x1.eff61cb926c5ep-2,
                            0x1.22100156b5d05p+12, -0x1.33c3366da76f2p+12,
                            0x1.8d94e9dad97e9p+9, 0x1.725851d6ce3d8p+13,
                            -0x1.cb01d465fc1e5p+13, -0x1.6dfc90611a7f1p+13,
                            0x1.a0d3c9cb465a1p+14, -0x1.1649974d86d5ep+12,
                            -0x1.284a3b2d44c7fp+14, 0x1.d97934642011cp+14,
                            -0x1.8c9583b0822a0p+13, -0x1.468f0def3f147p+15,
                            0x1.631c3028c9906p+15, 0x1.388320e76051dp+14,
                            -0x1.7296007052ad3p+15, 0x1.6a1fde4d7bb1cp+14,
                        },
                    },
                },
            };
            // clang-format on
        } // namespace anchor

        namespace taylor
        {
            // clang-format off
            alignas(64) const signed char Tiles[5000] = {
            -2,-2, 0,16, 0,17, 1,17, 1,14, 1,18, 2,19, 2,16, 2,14, 2,16, 2,21, 3,18, 3,15, 3,13, 3,15, 3,18, 4,17, 4,15, 4,13, 4,13, 4,15, 4,18, 4,20, 5,21, 5,19, 5,17, 5,16, 5,14, 5,12, 5,12, 5,14, 5,16, 5,18, 5,20, 5,22, 6,20, 6,19, 6,18, 6,16, 6,15, 6,13, 6,12, 6,10, 6,12, 6,13, 6,15, 6,16, 6,18, 6,19,-2,-2,
             7,15, 7,15, 8,16, 8,16, 1,18, 9,14, 9,17, 2,19, 2,17, 2,19,12,18, 3,19, 3,17, 3,16, 3,17, 3,20, 4,19, 4,17, 4,15, 4,15, 4,17, 4,19, 4,22, 5,21, 5,20, 5,18, 5,16, 5,15, 5,14, 5,14, 5,15, 5,17, 5,18, 5,20, 5,22, 6,21, 6,19, 6,18, 6,17, 6,15, 6,14, 6,13, 6,12, 6,13, 6,14, 6,15, 6,17, 6,18, 6,20,-2,-2,
            10,15,10,16, 8,16, 8,16, 9,17, 9,14, 9,17,11,16,11,18,12,18,12,16,12,16, 3,20, 3,19, 3,20,16,20, 4,20, 4,19, 4,18, 4,18, 4,19, 4,21,17,19,17,21, 5,21, 5,19, 5,18, 5,17, 5,16, 5,16, 5,17, 5,18, 5,19, 5,21, 6,22, 6,21, 6,20, 6,18, 6,17, 6,16, 6,15, 6,14, 6,14, 6,14, 6,15, 6,16, 6,17, 6,19, 6,20,-2,-2,
            13,15,14,16,14,17,15,17,15,14,15,17,11,16,11,12,11,16,12,16,12,13,12,13,12,17,16,17,16,16,16,17,16,19, 4,21, 4,20, 4,20, 4,21,17,16,17,17,17,19, 5,22, 5,20, 5,19, 5,18, 5,18, 5,18, 5,18, 5,19, 5,21, 5,23,28,21, 6,22, 6,20, 6,19, 6,18, 6,17, 6,16, 6,16, 6,15, 6,16, 6,16, 6,17, 6,18, 6,20, 6,21,-2,-2,
            13,16,14,16,14,17,15,17,15,14,15,18,11,18,11,16,11,18,12,18,12,16,12,16,12,18,16,15,16,13,16,15,16,18,17,20,17,17,17,15,17,14,17,14,17,16,17,18,17,20, 5,22, 5,21, 5,20, 5,19, 5,19, 5,20, 5,21, 5,22,28,20,28,20, 6,22, 6,21, 6,20, 6,19, 6,18, 6,18, 6,17, 6,17, 6,17, 6,18, 6,18, 6,19, 6,21, 6,22,-2,-2,
            18,16,19,16,19,17,20,18,20,16,20,19,21,18,21,16,21,16,21,19,12,19,12,19,16,17,16,15,16,13,16,15,16,18,17,19,17,17,17,14,17,12,17,12,17,14,17,17,17,20,27,19,27,19, 5,22, 5,21, 5,21, 5,22,28,20,28,19,28,18,28,18,28,18, 6,22, 6,21, 6,20, 6,20, 6,19, 6,18, 6,18, 6,19, 6,19, 6,20, 6,21, 6,22, 6,23,-2,-2,
            18,16,19,16,19,17,20,16,20,12,20,16,21,17,21,14,21,14,21,17,24,18,24,18,16,19,16,17,16,16,16,17,16,19,17,20,17,17,17,15,17,14,17,14,17,16,17,18,27,18,27,17,27,17,27,18,27,19,27,21,28,20,28,18,28,17,28,17,28,16,28,17,28,17, 6,22, 6,21, 6,21, 6,20, 6,20, 6,20, 6,20, 6,20, 6,21, 6,22, 6,23,32,20,-2,-2,
            22,15,22,15,23,18,23,18,20,16,20,19,21,19,21,16,21,16,21,19,24,15,24,15,24,18,16,19,16,19,16,20,31,21,31,21,17,19,17,17,17,16,17,16,17,17,17,19,27,16,27,15,27,15,27,16,27,18,27,19,28,19,28,17,28,16,28,15,28,15,28,15,28,16,28,17,28,19, 6,22, 6,22, 6,21, 6,21, 6,21, 6,22, 6,22,32,17,32,18,32,19,-2,-2,
            25,16,25,17,23,15,23,15,26,19,26,17,26,17,26,19,21,20,24,16,24,13,24,13,24,16,24,20,31,21,31,20,31,19,31,19,31,20,17,19,17,18,17,18,17,19,27,16,27,15,27,13,27,13,27,15,27,16,27,18,28,18,28,16,28,15,28,14,28,13,28,14,28,15,28,16,28,18,28,20,32,17,32,16, 6,23,32,14,32,14,32,15,32,16,32,17,32,18,-2,-2,
            25,17,25,17,23,18,23,19,26,17,26,14,26,14,26,18,24,20,24,18,24,15,24,15,24,18,24,21,31,19,31,18,31,17,31,17,31,18,31,20,17,21,17,21,27,18,27,16,27,14,27,11,27,11,27,14,27,16,27,18,28,17,28,16,28,14,28,12,28,11,28,12,28,14,28,16,28,17,32,17,32,16,32,15,32,14,32,13,32,13,32,14,32,15,32,16,32,17,-2,-2,
            29,16,29,16,30,17,30,17,26,19,26,17,26,17,33,18,33,20,24,20,24,18,24,18,24,20,31,20,31,18,31,16,31,14,31,14,31,16,31,18,31,21,27,20,27,18,27,16,27,15,27,13,27,13,27,15,27,17,27,18,28,17,28,16,28,14,28,12,28,11,28,12,28,14,28,16,28,17,32,17,32,15,32,14,32,13,32,11,32,11,32,13,32,14,32,15,-2,-2,-2,-2,
            29,17,29,17,30,14,30,14,30,19,33,18,33,15,33,15,33,18,36,19,36,18,36,18,36,19,31,20,31,17,31,15,31,12,31,12,31,15,31,18,31,20,27,21,27,19,27,17,27,16,27,15,27,15,27,16,27,18,27,19,28,18,28,16,28,15,28,14,28,13,28,14,28,15,28,16,28,18,32,16,32,15,32,14,32,12,32,10,32,10,32,12,32,14,32,15,-2,-2,-2,-2,
            34,18,34,18,30,18,35,16,35,18,33,18,33,15,33,15,33,18,36,17,36,15,36,15,36,17,36,20,31,18,31,16,31,14,31,14,31,16,31,18,31,21,39,22,27,20,27,19,27,18,27,17,27,17,27,18,27,19,27,21,28,19,28,17,28,16,28,15,28,15,28,15,28,16,28,17,28,19,32,17,32,15,32,14,32,13,32,11,32,11,32,13,32,14,32,15,-2,-2,-2,-2,
            34,15,34,15,35,16,35,12,35,16,33,20,33,18,33,18,33,20,36,16,36,13,36,13,36,16,36,19,31,20,31,18,31,17,31,17,31,18,31,20,31,23,39,20,39,20,27,21,27,20,27,19,27,19,27,20,27,21,28,21,28,20,28,18,28,17,28,17,28,16,28,17,28,18,28,19,28,20,32,17,32,16,32,15,32,14,32,13,32,13,32,14,32,15,32,16,-2,-2,-2,-2,
            34,18,34,18,35,18,35,16,35,19,38,20,38,19,38,19,38,20,36,17,36,15,36,15,36,17,36,20,31,21,31,20,31,19,31,19,31,20,39,18,39,17,39,18,39,18,39,20,39,21,27,21,27,22,43,22,43,22,43,22,28,21,28,20,28,19,28,18,28,18,28,18,28,19,28,20,28,22,32,18,32,17,32,16,32,15,32,14,32,14,32,15,32,16,-2,-2,-2,-2,-2,-2,
            37,20,37,17,37,17,37,20,38,21,38,18,38,16,38,16,38,18,36,19,36,18,36,18,36,20,40,21,40,21,40,22,31,22,31,22,39,18,39,17,39,15,39,15,39,17,39,18,39,20,39,22,43,21,43,21,43,20,43,20,43,20,43,21,28,21,28,20,28,20,28,20,28,21,28,22,32,20,32,19,32,18,32,17,32,16,32,16,32,16,32,16,32,17,-2,-2,-2,-2,-2,-2,
            37,18,37,14,37,14,37,18,38,20,38,17,38,13,38,13,38,17,38,20,36,21,36,21,40,20,40,19,40,19,40,20,40,21,39,19,39,17,39,15,39,13,39,13,39,15,39,17,39,19,39,21,43,20,43,19,43,19,43,18,43,19,43,19,43,20,43,21,28,22,45,22,45,22,45,22,32,21,32,20,32,19,32,18,32,18,32,17,32,17,32,18,32,18,-2,-2,-2,-2,-2,-2,
            37,20,37,17,37,17,37,21,38,21,38,18,38,16,38,16,38,18,38,21,40,21,40,20,40,18,40,17,40,17,40,18,40,20,39,18,39,16,39,14,39,12,39,12,39,14,39,16,39,19,39,21,43,19,43,18,43,17,43,17,43,17,43,18,43,19,43,20,43,22,45,21,45,20,45,20,45,21,45,21,32,20,32,20,32,19,32,19,32,19,32,19,-2,-2,-2,-2,-2,-2,-2,-2,
            41,21,41,19,41,17,41,19,41,22,38,21,38,19,38,19,38,21,42,22,40,20,40,18,40,16,40,14,40,14,40,16,40,18,39,19,39,17,39,15,39,13,39,13,39,15,39,17,39,19,39,21,43,17,43,16,43,15,43,15,43,15,43,16,43,18,43,19,43,21,45,19,45,19,45,19,45,19,45,20,45,20,45,21,32,20,32,20,32,20,32,21,-2,-2,-2,-2,-2,-2,-2,-2,
            41,20,41,17,41,14,41,17,41,20,42,19,42,17,42,16,42,17,42,19,40,20,40,17,40,15,40,12,40,12,40,15,40,18,40,20,39,18,39,17,39,15,39,16,39,17,39,18,39,20,43,18,43,16,43,15,43,14,43,13,43,14,43,15,43,17,43,18,43,20,45,18,45,18,45,17,45,18,45,18,45,19,45,20,45,21,45,22,32,22,32,22,-2,-2,-2,-2,-2,-2,-2,-2,
            41,20,41,17,41,14,41,17,41,20,42,18,42,15,42,13,42,15,42,18,42,21,40,18,40,16,40,14,40,14,40,16,40,18,40,21,39,20,39,18,39,18,39,18,39,19,39,20,39,22,43,17,43,16,43,14,43,12,43,11,43,12,43,14,43,16,43,18,45,18,45,17,45,16,45,16,45,16,45,17,45,18,45,19,45,20,45,21,45,23,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            41,22,41,19,41,18,41,19,41,22,42,18,42,15,42,13,42,15,42,18,42,21,40,20,40,18,40,17,40,17,40,18,40,20,40,22,39,21,39,20,39,20,39,20,39,21,39,22,43,19,43,17,43,16,43,14,43,12,43,11,43,12,43,14,43,16,43,18,45,17,45,16,45,15,45,15,45,15,45,16,45,17,45,18,45,19,45,20,45,22,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            44,23,44,20,44,19,44,19,44,20,42,20,42,17,42,16,42,17,42,20,46,19,46,20,40,20,40,19,40,19,40,20,40,22,47,22,47,22,39,23,39,22,39,22,39,23,43,21,43,20,43,18,43,17,43,15,43,14,43,13,43,14,43,15,43,17,43,18,45,16,45,15,45,14,45,13,45,14,45,15,45,16,45,17,45,18,45,20,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            44,21,44,18,44,16,44,16,44,18,44,21,42,20,42,19,46,18,46,17,46,17,46,18,46,20,40,22,40,22,40,23,47,21,47,20,47,20,47,20,47,20,47,21,47,23,43,22,43,20,43,19,43,18,43,16,43,15,43,15,43,15,43,16,43,18,45,16,45,15,45,14,45,12,45,11,45,12,45,14,45,15,45,16,45,18,45,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            44,20,44,17,44,13,44,13,44,17,44,20,46,21,46,18,46,16,46,14,46,14,46,16,46,18,46,21,47,22,47,21,47,19,47,18,47,18,47,18,47,19,47,20,47,21,47,23,43,21,43,20,43,19,43,18,43,17,43,17,43,17,43,18,43,19,45,16,45,15,45,13,45,11,45, 9,45,11,45,13,45,15,45,16,45,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            44,21,44,18,44,16,44,16,44,18,44,21,46,20,46,18,46,15,46,12,46,12,46,15,46,18,46,20,47,21,47,19,47,18,47,17,47,16,47,16,47,17,47,18,47,20,47,21,43,23,43,21,43,20,43,19,43,19,43,18,43,19,43,19,45,18,45,16,45,15,45,14,45,12,45,11,45,12,45,14,45,15,45,16,45,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            44,23,44,21,44,19,44,19,44,21,44,23,46,21,46,18,46,16,46,14,46,14,46,16,46,19,46,21,47,20,47,18,47,17,47,15,47,14,47,14,47,15,47,17,47,19,47,21,47,23,43,23,43,22,43,21,43,20,43,20,43,21,43,21,45,18,45,17,45,16,45,15,45,14,45,13,45,14,45,15,45,16,45,17,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,23,48,21,48,19,48,18,48,19,48,21,48,23,46,20,46,18,46,17,46,17,46,18,46,20,46,22,47,20,47,18,47,16,47,14,47,12,47,12,47,14,47,16,47,18,47,20,47,22,49,21,49,21,49,20,49,20,43,22,43,22,45,20,45,19,45,18,45,17,45,16,45,15,45,15,45,15,45,16,45,17,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,21,48,19,48,17,48,15,48,17,48,19,48,21,46,22,46,20,46,19,46,19,46,20,46,22,47,22,47,20,47,18,47,16,47,14,47,12,47,12,47,14,47,16,47,18,47,20,47,22,49,20,49,19,49,19,49,19,49,19,49,19,49,20,45,20,45,19,45,18,45,17,45,16,45,16,45,16,45,17,45,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,20,48,17,48,15,48,12,48,15,48,17,48,20,48,23,46,23,46,22,46,22,46,23,50,23,47,22,47,20,47,19,47,17,47,15,47,14,47,14,47,15,47,17,47,19,47,21,47,23,49,19,49,18,49,17,49,17,49,17,49,18,49,19,49,20,45,20,45,19,45,18,45,18,45,18,45,18,45,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,20,48,17,48,15,48,12,48,15,48,17,48,20,48,23,50,22,50,21,50,20,50,20,50,20,50,21,47,21,47,20,47,18,47,17,47,16,47,16,47,17,47,18,47,20,47,22,49,19,49,18,49,17,49,16,49,16,49,16,49,17,49,18,49,19,49,20,45,20,45,20,45,19,45,19,45,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,21,48,19,48,17,48,15,48,17,48,19,48,21,50,22,50,21,50,19,50,18,50,18,50,18,50,19,50,21,47,21,47,20,47,19,47,18,47,18,47,19,47,20,47,21,49,19,49,18,49,17,49,15,49,15,49,14,49,15,49,15,49,17,49,18,49,19,49,21,45,21,45,21,45,21,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            48,23,48,21,48,19,48,19,48,19,48,21,48,23,50,21,50,19,50,18,50,16,50,16,50,16,50,18,50,19,50,21,47,22,47,21,47,20,47,20,47,21,47,22,47,23,49,18,49,17,49,16,49,14,49,13,49,12,49,13,49,14,49,16,49,17,49,19,49,20,49,22,45,22,45,22,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,23,51,21,48,22,48,22,48,22,51,23,50,22,50,20,50,18,50,16,50,15,50,14,50,15,50,16,50,18,50,20,50,22,47,23,47,22,47,22,47,23,49,21,49,20,49,18,49,17,49,15,49,13,49,12,49,10,49,12,49,13,49,15,49,17,49,18,49,20,49,21,49,23,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,20,51,19,51,18,51,18,51,19,51,20,51,22,50,19,50,17,50,15,50,13,50,11,50,13,50,15,50,17,50,19,50,21,52,23,52,22,52,23,49,23,49,21,49,20,49,18,49,17,49,15,49,13,49,12,49,10,49,12,49,13,49,15,49,17,49,18,49,20,49,21,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,19,51,17,51,16,51,16,51,17,51,19,51,21,50,19,50,17,50,15,50,13,50,11,50,13,50,15,50,17,50,19,50,21,52,21,52,21,52,21,52,21,49,21,49,20,49,18,49,17,49,16,49,14,49,13,49,12,49,13,49,14,49,16,49,17,49,19,49,20,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,17,51,15,51,13,51,13,51,15,51,17,51,20,50,20,50,18,50,16,50,15,50,14,50,15,50,16,50,18,50,20,50,22,52,19,52,19,52,19,52,20,52,21,49,21,49,19,49,18,49,17,49,15,49,15,49,14,49,15,49,16,49,17,49,18,49,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,16,51,14,51,11,51,11,51,14,51,16,51,19,51,22,50,19,50,18,50,17,50,16,50,17,50,18,50,19,50,21,52,18,52,18,52,17,52,18,52,18,52,19,52,20,49,20,49,19,49,18,49,17,49,16,49,16,49,16,49,17,49,18,49,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,17,51,15,51,13,51,13,51,15,51,17,51,20,51,22,50,21,50,20,50,19,50,18,50,19,50,20,50,21,52,18,52,17,52,16,52,16,52,16,52,17,52,18,52,19,52,20,49,20,49,19,49,18,49,18,49,17,49,18,49,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,19,51,17,51,16,51,16,51,17,51,19,51,21,53,22,50,23,50,22,50,21,50,20,50,21,50,22,52,18,52,17,52,16,52,15,52,14,52,15,52,16,52,17,52,18,52,20,52,21,49,20,49,20,49,19,49,19,49,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,21,51,19,51,18,51,18,51,19,51,21,51,23,53,21,53,21,53,22,50,23,50,23,50,23,52,19,52,17,52,16,52,14,52,13,52,12,52,13,52,14,52,16,52,17,52,19,52,21,49,22,49,21,49,21,49,21,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            51,23,51,22,51,21,51,21,51,22,53,18,53,18,53,19,53,20,53,21,53,22,52,22,52,20,52,18,52,17,52,15,52,14,52,12,52,10,52,12,52,14,52,15,52,17,52,19,52,20,52,22,49,23,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,22,53,21,53,19,53,18,53,17,53,16,53,16,53,17,53,18,53,19,53,21,53,23,52,20,52,18,52,17,52,15,52,14,52,12,52,10,52,12,52,14,52,15,52,17,52,19,52,20,52,22,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,21,53,20,53,18,53,17,53,15,53,15,53,15,53,15,53,17,53,18,53,20,53,21,52,20,52,19,52,17,52,16,52,14,52,13,52,12,52,13,52,14,52,16,52,17,52,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,21,53,19,53,17,53,15,53,14,53,13,53,13,53,14,53,15,53,17,53,19,53,21,52,21,52,20,52,18,52,17,52,16,52,15,52,14,52,15,52,16,52,17,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,20,53,18,53,16,53,15,53,13,53,11,53,11,53,13,53,15,53,16,53,18,53,20,53,22,52,21,52,19,52,18,52,17,52,16,52,16,52,16,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,20,53,18,53,16,53,15,53,13,53,11,53,11,53,13,53,15,53,16,53,18,53,20,53,22,52,22,52,20,52,19,52,19,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,21,53,19,53,17,53,16,53,14,53,13,53,13,53,14,53,16,53,17,53,19,53,21,53,23,52,23,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            53,22,53,20,53,18,53,17,53,16,53,15,53,15,53,16,53,17,53,18,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            };

            alignas(64) const double TaylorCoeffs[2 * 24 * 54] = {
            0x0.0p+0, 0x1.249249249248bp-2, 0x1.7d3335cdc4259p-1, 0x0.0p+0, 0x0.0p+0, 0x1.67e6cb499a33fp-1, -0x1.165efbb8bca93p-1, 0x0.0p+0, 0x0.0p+0, -0x1.75d296c66d5a5p-2, 0x1.c1ef8338c7cb0p-3, 0x0.0p+0, 0x0.0p+0, 0x1.ef4615133f783p-4, -0x1.f99386606c042p-5, 0x0.0p+0, 0x0.0p+0, -0x1.e37bed3fb983ap-6, 0x1.b481d232bf3a9p-7, 0x0.0p+0, 0x0.0p+0, 0x1.76558b7f680c2p-8, -0x1.326cb5791f60cp-9, 0x0.0p+0, 0x0.0p+0, -0x1.e0d03aceb2119p-11, 0x1.6ac641b262057p-12, 0x0.0p+0, 0x0.0p+0, 0x1.07fe2f4af19f7p-13, -0x1.737f9c2d583fep-15, 0x0.0p+0, 0x0.0p+0, -0x1.fa935948836cbp-17, 0x1.4f50838398015p-18, 0x0.0p+0, 0x0.0p+0, 0x1.afb155b5c8102p-20, -0x1.0ea5f58e00823p-21, 0x0.0p+0, 0x0.0p+0, -0x1.4af8b5391433ep-23, 0x1.8b3640f9efe06p-25, 0x0.0p+0, 0x0.0p+0, 0x1.cd523f9496f64p-27, -0x1.07764d330a5aap-28, 0x0.0p+0,
            0x0.0p+0, 0x1.4924924924939p-1, 0x1.192e96b2eb224p-1, 0x0.0p+0, 0x0.0p+0, 0x1.b06bbe44935bdp-2, -0x1.1c610e82e53f4p-2, 0x0.0p+0, 0x0.0p+0, -0x1.4ccec07d86128p-3, 0x1.62cf7cb518df8p-4, 0x0.0p+0, 0x0.0p+0, 0x1.5e04aab4d8598p-5, -0x1.4312a6426f1b7p-6, 0x0.0p+0, 0x0.0p+0, -0x1.19574c8c7ad92p-7, 0x1.d1489b4f21860p-9, 0x0.0p+0, 0x0.0p+0, 0x1.6f393f44a82e8p-10, -0x1.15cc47682a15cp-11, 0x0.0p+0, 0x0.0p+0, -0x1.94438dc2201e6p-13, 0x1.1bc4f198250a0p-14, 0x0.0p+0, 0x0.0p+0, 0x1.814b8f996a756p-16, -0x1.fb141d42bfd64p-18, 0x0.0p+0, 0x0.0p+0, -0x1.440dc27af1951p-19, 0x1.92eb1ebfcd9ccp-21, 0x0.0p+0, 0x0.0p+0, 0x1.e8179cb601029p-23, -0x1.206c194348e83p-24, 0x0.0p+0, 0x0.0p+0, -0x1.4cf4b30617760p-26, 0x1.77dbbfdf1e28fp-28, 0x0.0p+0, 0x0.0p+0, 0x1.9f507951c549bp-30, -0x1.c19ddccdadf4cp-32, 0x0.0p+0,
            0x0.0p+0, 0x1.36db6db6db718p+0, 0x1.806df05850832p-2, 0x0.0p+0, 0x0.0p+0, 0x1.bbb118bd7d6c1p-3, -0x1.cc2e3c275a8b4p-4, 0x0.0p+0, 0x0.0p+0, -0x1.b61fafd1f8c7ep-5, 0x1.845a90869941ep-6, 0x0.0p+0, 0x0.0p+0, 0x1.43bd76f916690p-7, -0x1.ff759c352b3c4p-9, 0x0.0p+0, 0x0.0p+0, -0x1.81164ff51a6e4p-10, 0x1.15a8103bc8e42p-11, 0x0.0p+0, 0x0.0p+0, 0x1.80e81b19d00b2p-13, -0x1.014bb509f99fep-14, 0x0.0p+0, 0x0.0p+0, -0x1.4ca45c51ee68fp-16, 0x1.a0d71e216af09p-18, 0x0.0p+0, 0x0.0p+0, 0x1.fb537b26dbfc1p-20, -0x1.2c60b903cfd8bp-21, 0x0.0p+0, 0x0.0p+0, -0x1.5a9e106b3f76dp-23, 0x1.864f31343926ap-25, 0x0.0p+0, 0x0.0p+0, 0x1.ad6c2535b799cp-27, -0x1.ce21a8fac58d6p-29, 0x0.0p+0, 0x0.0p+0, -0x1.e6f5c5a08336cp-31, 0x1.f6e30d7678aedp-33, 0x0.0p+0, 0x0.0p+0, 0x1.fd68095012b86p-35, -0x1.fa8a08d21ba02p-37, 0x0.0p+0,
            0x0.0p+0, 0x1.edb6db6db6dbfp+0, 0x1.0d8f964c3c1c3p-2, 0x0.0p+0, 0x0.0p+0, 0x1.cee6d3aa69b5ep-4, -0x1.73025d62db6b4p-5, 0x0.0p+0, 0x0.0p+0, -0x1.186175115f936p-6, 0x1.92918d6cfd05cp-8, 0x0.0p+0, 0x0.0p+0, 0x1.141cc84d9cd56p-9, -0x1.6b83454edead0p-11, 0x0.0p+0, 0x0.0p+0, -0x1.cd0434615f8abp-13, 0x1.1a79581b45caap-14, 0x0.0p+0, 0x0.0p+0, 0x1.4f57f6d7c06e1p-16, -0x1.82878e781cd6dp-18, 0x0.0p+0, 0x0.0p+0, -0x1.b166560daa40dp-20, 0x1.d984629c72014p-22, 0x0.0p+0, 0x0.0p+0, 0x1.f8da916d978f8p-24, -0x1.06fa4a5c48e65p-25, 0x0.0p+0, 0x0.0p+0, -0x1.0c04f46dbf84ep-27, 0x1.0b8216bc748c3p-29, 0x0.0p+0, 0x0.0p+0, 0x1.05b9fc4e932b0p-31, -0x1.f677531b6d46dp-34, 0x0.0p+0, 0x0.0p+0, -0x1.d9966501cce19p-36, 0x1.b69b5b9f85d19p-38, 0x0.0p+0, 0x0.0p+0, 0x1.8f6a531e16516p-40, -0x1.65dd3a5b3ea10p-42, 0x0.0p+0,
            0x0.0p+0, 0x1.5b6db6db6db56p+1, 0x1.914c2ad1ed440p-3, 0x0.0p+0, 0x0.0p+0, 0x1.08e57f5b66e68p-4, -0x1.4e5c821cc3d36p-6, 0x0.0p+0, 0x0.0p+0, -0x1.956c4b604655dp-8, 0x1.da05e56dc564ap-10, 0x0.0p+0, 0x0.0p+0, 0x1.0c0be91284ec1p-11, -0x1.25fcffede4340p-13, 0x0.0p+0, 0x0.0p+0, -0x1.3964537d28d33p-15, 0x1.4551663f46968p-17, 0x0.0p+0, 0x0.0p+0, 0x1.49635edb34113p-19, -0x1.45c64d54e5dcap-21, 0x0.0p+0, 0x0.0p+0, -0x1.3b21efa78be84p-23, 0x1.2a7e105de1b6bp-25, 0x0.0p+0, 0x0.0p+0, 0x1.15215186a44d6p-27, -0x1.f8de2959e0614p-30, 0x0.0p+0, 0x0.0p+0, -0x1.c390d1f6e10e9p-32, 0x1.8ce57ee4b3e85p-34, 0x0.0p+0, 0x0.0p+0, 0x1.570bbbf7e21dcp-36, -0x1.23c2b9c0e1dc6p-38, 0x0.0p+0, 0x0.0p+0, -0x1.e8a42ee491526p-41, 0x1.931bf39cf3283p-43, 0x0.0p+0, 0x0.0p+0, 0x1.47c6a45f72910p-45, -0x1.06d24f11b8171p-47, 0x0.0p+0,
            0x0.0p+0, 0x1.09249249248eep+2, 0x1.0f6716ed4a123p-3, 0x0.0p+0, 0x0.0p+0, 0x1.f13e2cbe1c2f4p-6, -0x1.bcdc8a7f6fa0fp-8, 0x0.0p+0, 0x0.0p+0, -0x1.854381a8edefep-10, 0x1.4d906accee48ep-12, 0x0.0p+0, 0x0.0p+0, 0x1.183be79a4b141p-14, -0x1.ce194ff9bc9d1p-17, 0x0.0p+0, 0x0.0p+0, -0x1.763e23aac905ep-19, 0x1.29f619da44021p-21, 0x0.0p+0, 0x0.0p+0, 0x1.d2c5d5935fb4ep-24, -0x1.67ef236f344d3p-26, 0x0.0p+0, 0x0.0p+0, -0x1.116858fd405eep-28, 0x1.9965acda2d5a0p-31, 0x0.0p+0, 0x0.0p+0, 0x1.2e44d33993e1ap-33, -0x1.b86123a5bfa0bp-36, 0x0.0p+0, 0x0.0p+0, -0x1.3ca6cb2fab085p-38, 0x1.c1acfbacf0decp-41, 0x0.0p+0, 0x0.0p+0, 0x1.3b6bd3e35ebd3p-43, -0x1.b54bbfc4243dap-46, 0x0.0p+0, 0x0.0p+0, -0x1.2bab67f894272p-48, 0x1.962805bf730e7p-51, 0x0.0p+0, 0x0.0p+0, 0x1.10451ee506fe7p-53, -0x1.69332159a4115p-56, 0x0.0p+0,
            0x0.0p+0, 0x1.84924924924b2p+2, 0x1.77a80f66f1becp-4, 0x0.0p+0, 0x0.0p+0, 0x1.e26260bb68034p-7, -0x1.31f5e43b7e28bp-9, 0x0.0p+0, 0x0.0p+0, -0x1.7f926e2ff28b6p-12, 0x1.db6aa5bec28a3p-15, 0x0.0p+0, 0x0.0p+0, 0x1.23639012b0080p-17, -0x1.61628e62afd35p-20, 0x0.0p+0, 0x0.0p+0, -0x1.a822e61b881ebp-23, 0x1.f7eed76ef5ef8p-26, 0x0.0p+0, 0x0.0p+0, 0x1.2871371fbd359p-28, -0x1.597243f839b9cp-31, 0x0.0p+0, 0x0.0p+0, -0x1.8ed09394bd1c6p-34, 0x1.c842bc9778f7cp-37, 0x0.0p+0, 0x0.0p+0, 0x1.02af7862cb501p-39, -0x1.22cd92f742f43p-42, 0x0.0p+0, 0x0.0p+0, -0x1.44273a7291a29p-45, 0x1.6659d97035a9fp-48, 0x0.0p+0, 0x0.0p+0, 0x1.88f69d0876a95p-51, -0x1.ab84e362095cep-54, 0x0.0p+0, 0x0.0p+0, -0x1.cd85200255816p-57, 0x1.ee741b17a3e14p-60, 0x0.0p+0, 0x0.0p+0, 0x1.06e7360ba5631p-62, -0x1.158a0667c9eb4p-65, 0x0.0p+0,
            0x1.b6db6db6db6dcp-3, 0x1.2492492492469p-3, 0x1.a60e11f1d7647p-1, 0x1.796bab313a581p-3, -0x1.33d81b7f1e00ep-2, 0x1.a0b4804718c89p-1, -0x1.498ad81a10a3ap-1, -0x1.435208a4eca6dp-2, 0x1.0c978c341a99ap-2, -0x1.bea5ee5aa6bc9p-2, 0x1.0cdc780e01e17p-2, 0x1.7ca91567c3708p-3, -0x1.de6c92032ad04p-4, 0x1.25f6a7d04fda2p-3, -0x1.284fdafe9c6abp-4, -0x1.10fc80d07f78dp-4, 0x1.1f6239f09e99cp-5, -0x1.1657b8425e3adp-5, 0x1.eb47461b9bac6p-7, 0x1.1a4807466bf45p-6, -0x1.04ec7b7544651p-7, 0x1.99df08d89953fp-8, -0x1.44df29bb17989p-9, -0x1.c8f64ebefb326p-9, 0x1.7d1a1326ea8f3p-10, -0x1.eb356352fe913p-11, 0x1.63540b9a9c7f7p-12, 0x1.30097109931a0p-11, -0x1.d1bf875dd21ebp-13, 0x1.ed2423263495ap-14, -0x1.48ee4ed7e558fp-15, -0x1.578e872f775abp-14, 0x1.e96b1801850e1p-16, -0x1.a66db7716a6bfp-17, 0x1.05680edcd34a9p-18, 0x1.51647fd02a066p-17, -0x1.c320986e6f8b6p-19, 0x1.37f5f2aeb070ap-20, -0x1.66fc9084e5672p-22, -0x1.2502b9a96f28ap-20, 0x1.725e936720c08p-22, -0x1.8e10f0a379fbdp-24, 0x1.a8a5e8c530f5cp-26, 0x1.c82edca64e069p-24, -0x1.121b3256826ecp-25, 0x1.b2825a67f68dfp-28, -0x1.a83a03473c372p-30, -0x1.41c32f5c88520p-27,
            0x1.24924924924cfp-2, 0x1.b6db6db6db6acp-2, 0x1.3e49bce90df95p-1, 0x1.3ff809eec4337p-3, -0x1.c541a7480a191p-3, 0x1.0333b46799242p-1, -0x1.5da6ebdcbbd0cp-2, -0x1.a5f26bd8937bbp-3, 0x1.3ad0c2f0e7334p-3, -0x1.9b0882b604c5ap-3, 0x1.b1331edd7cabep-4, 0x1.9476cfee7ae9dp-4, -0x1.d00c65624e573p-5, 0x1.a0adbb9977a70p-5, -0x1.72288ae8be430p-6, -0x1.e6114cc53695ep-6, 0x1.d7bcf52222174p-7, -0x1.3231f636d9130p-7, 0x1.da4ed3f95ee41p-9, 0x1.acb916d805272p-8, -0x1.6fc6aa1d686dep-9, 0x1.591edcdeb39a5p-10, -0x1.d8707f2690e72p-12, -0x1.2ba6fa6b5b925p-10, 0x1.d2175639dca72p-12, -0x1.2fd5f73fd938bp-13, 0x1.6d8cecb5cd8d5p-15, 0x1.5b6b967dba40dp-13, -0x1.f2176bc8c5bb3p-15, 0x1.9727b2b7f565fp-17, -0x1.9a8c5c86e458ap-19, -0x1.585ec702b1fc5p-16, 0x1.cc67843251ef6p-18, -0x1.62f444c0894d6p-21, 0x1.b5b7943176a3dp-24, 0x1.2a35de082a6f3p-19, -0x1.770482d16ab45p-21, -0x1.184d0492ce5e6p-29, 0x1.71859336aff7dp-27, -0x1.ca94b53b616d2p-23, 0x1.110a170071209p-24, 0x1.a6c0dd0bd208ep-28, -0x1.6935702b356dbp-29, 0x1.3d09cd4ac5d94p-26, -0x1.67626187e7a80p-28, -0x1.0d97be8d80b59p-30, 0x1.72098e488060fp-32, -0x1.8e22e1434820bp-30,
            0x1.24924924924a8p-2, 0x1.9249249249257p-1, 0x1.e82b57742f093p-2, 0x1.90965b317a95fp-4, -0x1.e65238973d2d3p-4, 0x1.4b1c3d9c26da1p-2, -0x1.828c4601f3b3dp-3, -0x1.8ce3be488935bp-4, 0x1.07936f029f46cp-4, -0x1.9268677829961p-4, 0x1.7d9c6424bdca3p-5, 0x1.30c496e8cecdbp-5, -0x1.3d60dead75082p-6, 0x1.4e5252bd184fdp-6, -0x1.113d9d7c232abp-7, -0x1.2fcaad79798adp-7, 0x1.0eeebfdd6f0b5p-8, -0x1.a39b79ab5e03bp-9, 0x1.3039efc7f886dp-10, 0x1.c6a6d5944dd64p-10, -0x1.6985cb7f84b5cp-11, 0x1.a206a5425fb21p-12, -0x1.10bfaf7ba9009p-13, -0x1.11fcb62dc9924p-12, 0x1.8d9fe2c01ce9dp-14, -0x1.525107def4c50p-15, 0x1.8ea3ff8bc9ba1p-17, 0x1.154968a7a3214p-15, -0x1.74d581fe5efb8p-17, 0x1.bd03d53e50c97p-19, -0x1.d3e47958d87f4p-21, -0x1.e4913af028f2cp-19, 0x1.3112cea43018ap-20, -0x1.c9f0d41b585a2p-23, 0x1.97396a0fe073cp-25, 0x1.74e4ffd4b7a15p-22, -0x1.bb391b41070c0p-24, 0x1.3656d1a326881p-27, -0x1.4aa5a321e07cfp-30, -0x1.0089625798249p-25, 0x1.219d3f35d18ebp-27, 0x1.9cd49017108a9p-35, -0x1.00f41bf666419p-33, 0x1.3f40b3045e795p-29, -0x1.5802a47398c5ep-31, -0x1.033114fa26257p-34, 0x1.8ec04de314e53p-36, -0x1.6ab67b5d1773ep-33,
            0x1.6db6db6db6d78p-2, 0x1.2492492492453p-3, 0x1.89c8cf2394ed7p-1, 0x1.2bf4aa672c2ecp-2, -0x1.dcd8f739ecaebp-2, 0x1.6617dfb4bdae1p-1, -0x1.017a1551b424bp-1, -0x1.e79d9496dafe0p-2, 0x1.8a111ab90f794p-2, -0x1.3850c5e2142c9p-2, 0x1.4999604104755p-3, 0x1.0f6e472694906p-2, -0x1.4b518165ba354p-3, 0x1.32eedd11bc416p-4, -0x1.f8ed1e617fd93p-6, -0x1.6ee4a9d06f559p-4, 0x1.7675ed72ab0f8p-5, -0x1.68d690c4cb9a8p-7, 0x1.a75c42955bfa6p-9, 0x1.64413f5107a26p-6, -0x1.3ea272e7cf3c7p-7, 0x1.44721c85a0940p-11, 0x1.15e201ebdf8fep-14, -0x1.0db16b6d74d7cp-8, 0x1.b25214e6b7deap-10, 0x1.40d4ebb463ce8p-13, -0x1.bcacbfdbfa459p-14, 0x1.4e226f4e2cfabp-11, -0x1.ecf5299d1c657p-13, -0x1.d7a2258b1edf6p-15, 0x1.b40e8789212c0p-16, -0x1.5db74077df2cap-14, 0x1.de69717aaa2f4p-16, 0x1.7026b4d429b39p-17, -0x1.224f9192e4575p-18, 0x1.3c31c20e973c7p-17, -0x1.949e0caf446f1p-19, -0x1.b1432ac703bd0p-20, 0x1.3486514addd6bp-21, -0x1.f60c240724ab7p-21, 0x1.2e6f14cf3178dp-22, 0x1.a5ca5e4b9eecfp-23, -0x1.15f932af1a658p-24, 0x1.623a655adf9dcp-24, -0x1.93c122f8a41e0p-26, -0x1.626e4175eb67ep-26, 0x1.b65940e3aa7ddp-28, -0x1.c03c9667f2f1dp-28,
            0x1.000000000001dp-1, 0x1.124924924922ap+0, 0x1.8127a6931003fp-2, 0x1.e2fcc1aede208p-4, -0x1.f9a58d8c64258p-4, 0x1.a2c358d39a0d2p-3, -0x1.86725e52e877cp-4, -0x1.67fc70909f30dp-4, 0x1.a43a426d00293p-5, -0x1.3eee02b2c0045p-5, 0x1.ca237a3b0640bp-7, 0x1.ad3176bcb0852p-6, -0x1.8c1ca7ec961a8p-7, 0x1.1c87c488662a5p-8, -0x1.1f14ebc3aa645p-10, -0x1.50be0413fee4fp-8, 0x1.0b095617f47b8p-9, -0x1.66f8485b5831fp-13, -0x1.c07f4051831b6p-16, 0x1.8e978075c0526p-11, -0x1.19cfb28c5c7f6p-12, -0x1.50217dec47294p-15, 0x1.8c8a9ccd727abp-16, -0x1.7b6345f847479p-14, 0x1.e823802499458p-16, 0x1.6de168a5f8825p-17, -0x1.2866fa10bf0bcp-18, 0x1.2d066ce6fe42dp-17, -0x1.64b526ffe11abp-19, -0x1.b8217084d4978p-20, 0x1.31ed938c66a47p-21, -0x1.96ca05368d4f4p-21, 0x1.bef1cc39e922bp-23, 0x1.93223d7a53e5cp-23, -0x1.fba0b013760fap-25, 0x1.d95c3b5d36781p-25, -0x1.e3324c6d0fbd2p-27, -0x1.3322d611ee782p-26, 0x1.66a867b35aaa2p-28, -0x1.daee4fa6c2a83p-29, 0x1.c089464b47d3ep-31, 0x1.9579593d28892p-30, -0x1.bcee7bbbd3ae0p-32, 0x1.9578539042f81p-33, -0x1.5c7cf1d0ee522p-35, -0x1.dae06b5eb4489p-34, 0x1.edd649af81aaap-36, -0x1.19514e0e34099p-37,
            0x1.fffffffffffeap-2, 0x1.9249249249230p+0, 0x1.2ca7de9093c77p-2, 0x1.2b543fb4c55dbp-4, -0x1.05dfd69c751c5p-4, 0x1.0f6b7f034a29ep-3, -0x1.b54db665c4b1bp-5, -0x1.3e76daf21cb7ep-5, 0x1.429bdef5f92a2p-6, -0x1.3f021b5df4ab8p-6, 0x1.a80212553d3d4p-8, 0x1.2179d2ec19f28p-7, -0x1.da2066b7b2c3bp-9, 0x1.0034f6e00c4bbp-9, -0x1.15ca17e8bd0b5p-11, -0x1.68a249e775a3dp-10, 0x1.01b56576f09dbp-11, -0x1.0430f8d853d16p-13, 0x1.75bcad9db0c55p-16, 0x1.5cd2380a7ba6fp-13, -0x1.c1d16d459627ap-15, 0x1.5fcda2be13ec7p-20, 0x1.6206ae63a271fp-20, -0x1.1591e1a79796ep-16, 0x1.490319161eab7p-18, 0x1.dd3f30c9e7325p-21, -0x1.a8b6b91aba8ffp-22, 0x1.779f2386f2196p-20, -0x1.9df79700aa3dap-22, -0x1.3f7c6cabdcbefp-23, 0x1.b30decd2f81ccp-25, -0x1.b932809de9226p-24, 0x1.c75b79d119008p-26, 0x1.1421134731d08p-26, -0x1.4bef60b3e05e3p-28, 0x1.c794dcb3722d3p-28, -0x1.ba23e68038210p-30, -0x1.7d80f0700ef3cp-30, 0x1.a5e9c56f88c86p-32, -0x1.a059033636cfcp-32, 0x1.7c5af094e5a7ep-34, 0x1.c3062829a7936p-34, -0x1.d3a19be9405e4p-36, 0x1.50deff7bca678p-36, -0x1.20da1b39da450p-38, -0x1.d77f8535a68f7p-38, 0x1.cf4831387c946p-40, -0x1.de6b0f5ba26ffp-41,
            0x1.24924924924b8p-1, 0x1.24924924924cdp-4, 0x1.5e87a7031b62ep-1, 0x1.dd73c741c5123p-2, -0x1.6e8079c7fef28p-1, 0x1.fda62e56e25b7p-2, -0x1.eb98115b7427fp-3, -0x1.662955e247dbep-1, 0x1.1219059ba56a6p-1, -0x1.b969930564f29p-5, -0x1.1d984e7561ce2p-5, 0x1.62590872c8ca9p-2, -0x1.91fcfcb140be6p-3, -0x1.cf39b52426c02p-5, 0x1.8672259b07dfbp-5, -0x1.993467e27ed9cp-4, 0x1.7af77473c623ap-5, 0x1.01ffec6b3ec8cp-5, -0x1.2648eba19ee15p-6, 0x1.41f4bb9b25d53p-6, -0x1.f7cc7c30dc607p-8, -0x1.2dc1dddf0fad2p-7, 0x1.1c621c216bf53p-8, -0x1.6aca0505e5ebfp-9, 0x1.dd94df7e1a088p-11, 0x1.f3880322b202ap-10, -0x1.9cc49c3dd7541p-11, 0x1.1a0cdb9c3743ep-12, -0x1.1d18b20931c43p-14, -0x1.4320af572cf7cp-12, 0x1.e1d1b4b6ce856p-14, -0x1.a602ba8cd6acfp-17, 0x1.99d2e29a64763p-23, 0x1.5787c27491a51p-15, -0x1.d601f32a3e92ep-17, -0x1.6ba8abe2eb81ap-20, 0x1.e7638c22eef6bp-21, -0x1.354fbacc25696p-18, 0x1.887e12d8e1a3ep-20, 0x1.d41c6974b386bp-22, -0x1.830bb033dfd30p-23, 0x1.e0fbf806ea0dep-22, -0x1.1d0181773dfa3p-23, -0x1.23a78895eb7c6p-24, 0x1.9b0a3a073ffb9p-26, -0x1.46f6d7dccebc7p-25, 0x1.6b5b21b793f2bp-27, 0x1.12b6f9d3570dcp-27,
            0x1.249249249246ep-1, 0x1.2492492492485p-2, 0x1.2d4514a651b2cp-1, 0x1.5c8b2f8a75555p-2, -0x1.e972ef762e87ap-2, 0x1.9cd0f3d8e7fabp-2, -0x1.99d06c71b1a5cp-3, -0x1.bc989b0935bc6p-2, 0x1.3fac96c4369eep-2, -0x1.0b3e6882f9a58p-4, 0x1.92da04419b70ap-10, 0x1.877035bf49cb3p-3, -0x1.a772214c9fc49p-4, -0x1.230f759b5230ep-6, 0x1.1e833a45a6457p-6, -0x1.9d8a97a4b609dp-5, 0x1.71a4bfddb7324p-6, 0x1.859c492170927p-7, -0x1.ba142fbf91273p-8, 0x1.3113965bb05c0p-7, -0x1.d363edf504496p-9, -0x1.bd2497916e414p-9, 0x1.999f26b202196p-10, -0x1.4cf59cca7cd4ap-10, 0x1.b86bf79691a9ep-12, 0x1.5e5e5c3edaf18p-11, -0x1.19a8595a74530p-12, 0x1.0c86f7db821e3p-13, -0x1.28c14d55eaf2bp-15, -0x1.acea0bfec6c01p-14, 0x1.371b1d15ac77dp-15, -0x1.1d382c853d592p-17, 0x1.a4b9a507619a6p-20, 0x1.afd53a142c1b9p-17, -0x1.1fc60e2b38c65p-18, 0x1.76e8e94e0aa5cp-24, 0x1.cacd09dbde375p-24, -0x1.715060c348a5cp-20, 0x1.c96b049e43dd2p-22, 0x1.3ef45030a3136p-24, -0x1.2983f0136e0b5p-25, 0x1.11e50f46b33a8p-23, -0x1.3da19c756bfa4p-25, -0x1.d7946d32210eap-27, 0x1.534d7971b38dfp-28, -0x1.65208ba024a5ep-27, 0x1.85a71ee4c9174p-29, 0x1.c821a1652be2bp-30,
            0x1.2492492492477p-1, 0x1.4924924924934p-1, 0x1.de4caf6fa23f8p-2, 0x1.b518b1fb91300p-3, -0x1.09a331826e5c9p-2, 0x1.22bbfa83a7db9p-2, -0x1.173651ea13203p-3, -0x1.abd490a130d59p-3, 0x1.1531a954488c8p-3, -0x1.a400618c9a251p-5, 0x1.aa85e73f1979bp-7, 0x1.35a26891615b1p-4, -0x1.34a47ff68c9b1p-5, -0x1.ee8ab9ceadd73p-14, 0x1.709dacfb34686p-9, -0x1.1833d8a27fb3fp-6, 0x1.d583bb1fcf9bdp-8, 0x1.2eb07c2f7591ep-9, -0x1.634f19808677fp-10, 0x1.6e3f911890731p-9, -0x1.0b95d719413dfp-10, -0x1.61923f05c0474p-11, 0x1.3ba750957535cp-12, -0x1.6f90d186b6dfdp-12, 0x1.db668093da7e9p-14, 0x1.03bb1b505a979p-13, -0x1.901c53a89f11cp-15, 0x1.21455f1ffdee3p-15, -0x1.49f91a7ee6adfp-17, -0x1.23630879ce42ep-16, 0x1.941a3706b97acp-18, -0x1.5daeca29007c3p-19, 0x1.51bdc30b31c1cp-21, 0x1.0c2e3124cfc18p-19, -0x1.5606a2e043657p-21, 0x1.1c505d769ef9fp-23, -0x1.6bb011900f005p-26, -0x1.a47428e927510p-23, 0x1.f3764b4064ee1p-25, -0x1.a5d2ce1b39907p-31, -0x1.6b1db7d54f80fp-30, 0x1.1f4349bca36f9p-26, -0x1.4083c3df0bfb3p-28, -0x1.b469742fd1b87p-31, 0x1.6bfab8903c4fcp-32, -0x1.5b5d84a1dd3f5p-30, 0x1.6e17c91719be0p-32, 0x1.048305faa7ce2p-33,
            0x1.6db6db6db6defp-1, 0x1.0924924924917p+1, 0x1.d5e928bbeb808p-3, 0x1.13ebaecb83495p-4, -0x1.8f014077020dap-5, 0x1.4e1f4e2f652d7p-4, -0x1.a4d661f5c9a01p-6, -0x1.954c21c86eeadp-6, 0x1.596af50a0d173p-7, -0x1.d359a78f9c80ap-8, 0x1.bba0543e9f387p-10, 0x1.05ff4cf87aa09p-8, -0x1.6bc0e743ab084p-10, 0x1.430909c56029ap-12, -0x1.605873b8760f0p-16, -0x1.d5944ede4673fp-12, 0x1.1c98b0fa3bc7fp-13, 0x1.091e52d193b91p-16, -0x1.6a5f7c4c8485fp-17, 0x1.45f75575c8765p-15, -0x1.622ff4ee19984p-17, -0x1.3bf90ee0c6dcdp-18, 0x1.ccf8586dd3a6fp-20, -0x1.6ddfec5da938ep-19, 0x1.67715ac535ccbp-21, 0x1.2e6ffe06f4583p-21, -0x1.702366e6bdd57p-23, 0x1.4f464f3b72533p-23, -0x1.27948577c102dp-25, -0x1.a68f4c07ee233p-25, 0x1.cdee5f31ad59ap-27, -0x1.e7ff1f31b38c1p-28, 0x1.7213104cdc01fp-30, 0x1.e416af2074ab1p-29, -0x1.e891cbeae496ep-31, 0x1.edd929f93b33bp-33, -0x1.fdd89288a9220p-36, -0x1.dc790bc76b2dbp-33, 0x1.c22479cee9282p-35, -0x1.31bbb16d5d3c7p-40, -0x1.230675b005d8bp-40, 0x1.9cbf9f3cab000p-37, -0x1.6fd7b36fc296cp-39, -0x1.21ff94f1827e3p-41, 0x1.89439ebb14ff1p-43, -0x1.3ef3b737be719p-41, 0x1.0d41ef055bae3p-43, 0x1.c8f836b5d79cbp-45,
            0x1.9249249249229p-1, 0x1.8000000000030p+1, 0x1.5b6d5a62f0f27p-3, 0x1.4cef98398c465p-5, -0x1.746cbb4dd3d4cp-6, 0x1.7e3eca51e804ep-5, -0x1.82a7040274342p-7, -0x1.2a899d825c5dcp-7, 0x1.9817f7176e07ep-9, -0x1.668e71b29cf67p-9, 0x1.2c973abf2252dp-11, 0x1.f72f41e820be8p-11, -0x1.1f3303d7d6970p-12, 0x1.b6b3853294feap-14, -0x1.effd2dc1d86b9p-17, -0x1.33ed2ae20b900p-14, 0x1.38d5df1b115a4p-16, -0x1.89a00229b2fb4p-21, -0x1.15a832afbf500p-21, 0x1.2edd45200045fp-18, -0x1.186cf6eb4ec87p-20, -0x1.2cf180a19dc97p-22, 0x1.ab0eadc7b0a88p-24, -0x1.f1b4dad17337cp-23, 0x1.a7a55f0095d5cp-25, 0x1.00423eb7773b9p-25, -0x1.16cd52e4c572ep-27, 0x1.59b1ca2ec9196p-27, -0x1.0dcf6e44e4685p-29, -0x1.1b3a73fb9e691p-29, 0x1.10dde97544148p-31, -0x1.90ed27c3d0226p-32, 0x1.19051776dabc8p-34, 0x1.f746342732487p-34, -0x1.bef4e9f6ad169p-36, 0x1.6d375442bfa8dp-37, -0x1.a7abb8537434bp-40, -0x1.7fdc6a22bf2e9p-38, 0x1.3fc31c6aeba15p-40, -0x1.8aee1e92ea79ap-43, 0x1.47d18bfd0996bp-47, 0x1.02eca17252282p-42, -0x1.9842802849cf4p-45, -0x1.d5de4aa20a00ap-49, 0x1.cb56c3e34e5c1p-50, -0x1.39bae7d6de6d0p-47, 0x1.d65a52669a5adp-50, 0x1.16dc840e44cbdp-51,
            0x1.b6db6db6db69cp-1, 0x1.2492492492453p-4, 0x1.e3c36f71e43d8p-2, 0x1.1ac0ac4f74c88p-1, -0x1.7642d931b2f04p-1, 0x1.d3a78c5220b30p-4, 0x1.4c5b2a3742b77p-3, -0x1.322022b0225cbp-1, 0x1.76e6b0b3fcd3ap-2, 0x1.0800892801beap-2, -0x1.d4aa22c5555c4p-3, 0x1.672f557a43771p-3, -0x1.01f26bb640344p-4, -0x1.40f5da7eccf0fp-3, 0x1.72dbbdd6291bfp-4, -0x1.8f45cc8738dc3p-7, -0x1.22452325a7a40p-8, 0x1.781eba9073dbcp-5, -0x1.55e0154a69b4cp-6, -0x1.b7f51d81104fbp-8, 0x1.3dff589459145p-8, -0x1.199334a282960p-7, 0x1.a5e200aff0d86p-9, 0x1.67fa34ad7375fp-9, -0x1.616ab04e6e40dp-10, 0x1.1e67a9d651360p-10, -0x1.5aeb6f82cd811p-12, -0x1.39667404555bep-11, 0x1.0058d2ac51d14p-12, -0x1.6720a6fcca1bep-14, 0x1.1411f80dbf7bap-16, 0x1.87aeee7ae9e64p-14, -0x1.19c066b160325p-15, 0x1.3b34fb2f1e841p-21, 0x1.9f45ac8849cf7p-20, -0x1.7fbb45a32e18ep-17, 0x1.f084ee3277dfap-19, 0x1.1f889e45a3d65p-20, -0x1.14de6e50ea310p-21, 0x1.31d44b8dadd1dp-20, -0x1.66f659dde49e5p-22, -0x1.c68e41290fc79p-23, 0x1.5213011ad9f6ep-24, -0x1.91259e58f6ef4p-24, 0x1.a9aedb58595e5p-26, 0x1.d41d97d31d3ffp-26, -0x1.322209be23c0fp-27, 0x1.aa860b1215ed5p-28,
            0x1.b6db6db6db717p-1, 0x1.2492492492492p-2, 0x1.c53f760788851p-2, 0x1.ade213d694ef5p-2, -0x1.09aca8a4939ebp-1, 0x1.3f09b8a6b3021p-3, 0x1.7e2abde5959dcp-5, -0x1.9eccb976d6722p-2, 0x1.efd9d9d6841b9p-3, 0x1.e65531eabc813p-4, -0x1.c314ac344e028p-4, 0x1.df7fd627365a1p-4, -0x1.7672fea6bd551p-5, -0x1.33617ce904ab7p-4, 0x1.5d28423e8198bp-5, -0x1.b16ed4efdc975p-7, 0x1.9439765e502bfp-10, 0x1.5b5c4a6bd6f24p-6, -0x1.365c9d11a0776p-7, -0x1.60f526e1e4664p-10, 0x1.68b1ce81b7c62p-10, -0x1.f917be4c53effp-9, 0x1.790545d15b4dep-10, 0x1.be4f647fdf02ep-11, -0x1.c2fab09e5e440p-12, 0x1.029d484306d3fp-11, -0x1.4495416f69b11p-13, -0x1.925e8c3bb42cap-13, 0x1.47c1c333433c1p-14, -0x1.6f47ebdf9f2e0p-15, 0x1.68df3273afd79p-17, 0x1.f04922a6df793p-16, -0x1.6109a023aebffp-17, 0x1.12cbb7724091bp-19, -0x1.351e978c19307p-23, -0x1.db4b4d24dfb6cp-19, 0x1.3036eff29e5a8p-20, 0x1.03fd35289a3aep-23, -0x1.7a5d58d616c06p-24, 0x1.73790b30cbe96p-22, -0x1.b188e913a50c6p-24, -0x1.62d750ef9ecb3p-25, 0x1.17b5b0a673224p-26, -0x1.e3fad373356dep-26, 0x1.024e69e312b61p-27, 0x1.8df1d9b6e9e0ep-28, -0x1.07438c3157649p-29, 0x1.0721841fb0abcp-29,
            0x1.db6db6db6db74p-1, 0x1.4924924924943p-1, 0x1.73f07cccd11c3p-2, 0x1.25ee13fd7be8ep-2, -0x1.38d5a97f9be2cp-2, 0x1.06c3d3d2d5746p-3, 0x1.81493226b2d2ap-9, -0x1.ada3dbe3d105ep-3, 0x1.ca8260128e660p-4, 0x1.60d80e17b6d1ep-5, -0x1.446279ad1a7abp-5, 0x1.90b42a5543e6ep-5, -0x1.1e9130def6624p-6, -0x1.9d1d4b3b1d7a0p-6, 0x1.b0dcd1022918fp-7, -0x1.37707c243749ep-8, 0x1.44f5bafc0ba23p-11, 0x1.8bbf2f5c0c5b8p-8, -0x1.44847b7ea139fp-9, -0x1.489c4dba9da25p-12, 0x1.5857289391009p-12, -0x1.e435018b1338ep-11, 0x1.4ac4c8a162166p-12, 0x1.9298a0c68721bp-13, -0x1.7bb139994649dp-14, 0x1.9db6e92b2e2e9p-14, -0x1.d5b22de7a1afbp-16, -0x1.3b16ddde96457p-15, 0x1.dcdc99c5b6f7fp-17, -0x1.d8754ffc94097p-18, 0x1.89fa3aa4256b6p-20, 0x1.4f30327659d98p-18, -0x1.ba824bf6b79ebp-20, 0x1.a1ae9106284d7p-23, 0x1.93e64eb01e083p-26, -0x1.1430e4bf7338bp-21, 0x1.4768f700519e0p-23, 0x1.0eec9ea80931bp-25, -0x1.146ade75bb3b3p-26, 0x1.71820ae5779dap-25, -0x1.8d5fab4f04f5ep-27, -0x1.ba43fd382ae5cp-28, 0x1.37de00953a570p-29, -0x1.96e55aeda7830p-29, 0x1.8b8db19d72d22p-31, 0x1.94864b2ec78e5p-31, -0x1.ec9fd4ce17c85p-33, 0x1.6ac9f18f77815p-33,
            0x1.db6db6db6db42p-1, 0x1.2492492492494p+0, 0x1.351c2db457e33p-2, 0x1.670f909413b30p-3, -0x1.4769a86ff4d06p-3, 0x1.ce0995629c749p-4, -0x1.91619643061f8p-6, -0x1.8ec8cf718917dp-4, 0x1.85aa88e79f596p-5, 0x1.e971d90555a7ep-9, -0x1.f8ba48d4bd66ep-8, 0x1.43d4b01d7a87bp-6, -0x1.d3504a60453b0p-8, -0x1.5c46bf4f44b20p-8, 0x1.686bccd62ba2dp-9, -0x1.23e69441e77a8p-9, 0x1.306652e69a6b2p-11, 0x1.3d8a06cd1fc5bp-10, -0x1.f344caacfaa6cp-12, 0x1.cc00efe1ae078p-14, -0x1.877bfdca0ab72p-19, -0x1.6660c1f53ade7p-13, 0x1.dbb6d02dfcb6ep-15, 0x1.5f4259fecc443p-17, -0x1.ccf4a42f91bc5p-18, 0x1.25ed0879958adp-16, -0x1.5300e983fc45ap-18, -0x1.a67a9a3a77834p-19, 0x1.48d13b461163fp-20, -0x1.6c61af776cb6bp-20, 0x1.6a6d645c783cfp-22, 0x1.cd70c74f344c0p-22, -0x1.2bda55d958757p-23, 0x1.47ea970f566a9p-24, -0x1.0309fe1d361a8p-26, -0x1.6e38302d6caa4p-25, 0x1.a8113bfb3dcdep-27, -0x1.3bc3838eb02c4p-29, 0x1.e4bca3c4173d7p-34, 0x1.d43c83bd52710p-29, -0x1.eed30da4d01a3p-31, -0x1.cb3f796e19375p-34, 0x1.0a9febe18e0fdp-34, -0x1.f5a689a6b0674p-33, 0x1.e87452d74d9a5p-35, 0x1.9dc781c3c4325p-36, -0x1.10dc84a976abep-37, 0x1.c8d214defe61ep-37,
            0x1.1249249249239p+0, 0x1.249249249248ap-3, 0x1.53f8627803aa8p-2, 0x1.04ff73249c642p-1, -0x1.21aedae3d73d7p-1, -0x1.e1e3f57455a90p-5, 0x1.102c850e2373ep-2, -0x1.76b0b24d1d875p-2, 0x1.384c73db2a47ap-3, 0x1.19df410b9ef20p-2, -0x1.8f359541e364dp-3, 0x1.9304d98814008p-6, 0x1.886441a8b1598p-6, -0x1.c08d511080114p-4, 0x1.a37d6a6bf722ep-5, 0x1.e7c37d80a7d4fp-6, -0x1.5d08a26299545p-6, 0x1.491247a930b0fp-6, -0x1.a20141318c2efp-8, -0x1.8b3d34f87e798p-7, 0x1.80aec112a5b3bp-8, -0x1.6c1d856dc2cfbp-10, -0x1.6e8cd478676f0p-16, 0x1.4d38c56942de7p-9, -0x1.0557df3aa2da2p-10, -0x1.fbcd0847a0c30p-13, 0x1.70e4b9fc301c5p-13, -0x1.7611050676b97p-12, 0x1.e8ceea8b65b10p-14, 0x1.82ab27e2347bfp-14, -0x1.58a6289ca6311p-15, 0x1.212f6d592f8d2p-15, -0x1.2e6b854f3354fp-17, -0x1.13b7d651a476ep-16, 0x1.95f5d3de5155fp-18, -0x1.055b5942f45f2p-19, 0x1.1ca378a7c8bc5p-22, 0x1.16cc3db511363p-19, -0x1.68184ac41b47fp-21, -0x1.2dc7ecc150cfcp-25, 0x1.958fc535c48edp-25, -0x1.b762b0f749abap-23, 0x1.fb974f3da67ecp-25, 0x1.a96286e3251dap-26, -0x1.5e1ac5d75b51fp-27, 0x1.158d2193e1721p-26, -0x1.1e89647637041p-28, -0x1.fb7511923a483p-29,
            0x1.36db6db6db6f4p+0, 0x1.b6db6db6db6d5p-2, 0x1.1da5b80df1118p-2, 0x1.80b3d5a759d3ap-2, -0x1.6bf85100210b2p-2, -0x1.7a862dadd5ea4p-6, 0x1.245ab6bb1e190p-3, -0x1.8ff9e38434935p-3, 0x1.0cb17748e4c37p-4, 0x1.0fcdd7040c2cfp-3, -0x1.530087935aa1ep-4, 0x1.8b3d2df573e02p-9, 0x1.da62bb4b8cc6ep-7, -0x1.4aa7d1cb67ecdp-5, 0x1.0585eb4d5a164p-6, 0x1.b2a41fbf68fabp-7, -0x1.07c84d0f3fc20p-7, 0x1.461a7ad2a25dep-8, -0x1.09b82b23d10d5p-10, -0x1.03ca9b80f38c4p-8, 0x1.b9ae30c2e3543p-10, 0x1.fb1819e9b06f7p-15, -0x1.a61ac8b32ebc1p-13, 0x1.4c40b00db5664p-11, -0x1.bdb16c0d901cdp-13, -0x1.1b939056e5697p-13, 0x1.189716f7d9e6ep-14, -0x1.088d5dd4c98fep-14, 0x1.0d19fb6e22857p-16, 0x1.d8b3e3a4042a6p-16, -0x1.64277f66dab41p-17, 0x1.a594c6dcfae35p-19, -0x1.02643a7d7ed6fp-22, -0x1.eb10c5219e6d9p-19, 0x1.39585c8003af7p-20, 0x1.791c32b76ae06p-23, -0x1.164bfabd46184p-23, 0x1.740c15462ad94p-22, -0x1.9b18918ed317fp-24, -0x1.0215ce1c17f30p-24, 0x1.8df7c46827730p-26, -0x1.a47841de66ab5p-26, 0x1.885fd5f984c19p-28, 0x1.12787f7a6a8afp-27, -0x1.5d27c1fea8c1dp-29, 0x1.42b8fe004dfebp-30, -0x1.a82c862b43236p-33, -0x1.a01115c2f5a20p-31,
            0x1.36db6db6db6c1p+0, 0x1.9249249249273p+0, 0x1.d4b943bf3909fp-3, 0x1.25760c385787ep-3, -0x1.b00b8e87ca2fdp-4, 0x1.f463fab4ce42bp-5, -0x1.3a83a025a8114p-8, -0x1.a79afcc901a99p-5, 0x1.4837c4c0614a6p-6, 0x1.94174178e8a30p-8, -0x1.4257ec983cdc7p-8, 0x1.a189fc74d9e67p-8, -0x1.a883babc2c4f2p-10, -0x1.43a675376d23dp-9, 0x1.068e73dd154e0p-10, -0x1.2167c5b78ef02p-12, -0x1.24886fe4bbfb4p-17, 0x1.70a5de7a77b77p-12, -0x1.cbd24474eb16ap-14, -0x1.1fbf05e18dcfap-15, 0x1.47da6a988f112p-16, -0x1.00954dedbd295p-15, 0x1.fbd1a0f298511p-18, 0x1.115ebcf25a3d2p-17, -0x1.84afa66fe4bb2p-19, 0x1.ada5e0fc58e0ap-20, -0x1.12c65eb103b22p-22, -0x1.efb26db424f96p-21, 0x1.22096e6a2fb49p-22, -0x1.987438fd43213p-27, -0x1.bd9497097f1ebp-27, 0x1.3aed2737f036fp-24, -0x1.3ebc755023380p-26, -0x1.04159b6686212p-27, 0x1.957996ba1335fp-29, -0x1.2c7687a8a52d7p-28, 0x1.05eece8429ea2p-30, 0x1.0a8e11759b6e9p-30, -0x1.3b9e0dcf69dd1p-32, 0x1.9eba762423e65p-33, -0x1.1d2821995bf2bp-35, -0x1.5a203a2b21724p-34, 0x1.64ace3c00d6cap-36, -0x1.282cc18af182bp-38, 0x1.fc6ca81d90bafp-44, 0x1.5c2957e08ffefp-38, -0x1.435bdf9eb649dp-40, -0x1.9495f179dfabdp-43,
            0x1.492492492495dp+0, 0x1.2492492492478p-3, 0x1.cae61c67e9a7dp-3, 0x1.eea97beb1936dp-2, -0x1.c0ae4a426e618p-2, -0x1.6c2b6f90e4894p-3, 0x1.4169ede03fa60p-2, -0x1.88e9532d7f53ep-3, 0x1.3a6dbf0c61d4dp-8, 0x1.032b6f223fbdcp-2, -0x1.22b4c9442d149p-3, -0x1.12edb5bbc50b1p-4, 0x1.1370f0501d04cp-4, -0x1.e01cb5367aa45p-5, 0x1.014e83de22635p-6, 0x1.6ad0c84ff21efp-5, -0x1.7bb1716d29b5bp-6, -0x1.5afe60b3ff2ffp-13, 0x1.ccf7416e51f2dp-9, -0x1.4df3b002ac4a0p-7, 0x1.f4e61e91f4381p-9, 0x1.73d84767a3d43p-9, -0x1.a4be43122ea0ap-10, 0x1.3a773c4a4b0b7p-10, -0x1.2e7c79c58ce81p-12, -0x1.8b99494f3e0d6p-11, 0x1.4679d59a4c9cfp-12, -0x1.073216a00e4e9p-15, -0x1.2c3ae15aabfc4p-16, 0x1.e4384f501a0afp-14, -0x1.45c773b1ce895p-15, -0x1.126d73202ae5dp-16, 0x1.24ff933b8b8a1p-17, -0x1.8d97308121e72p-17, 0x1.b2cd85466582fp-19, 0x1.fd17d0eeed980p-19, -0x1.863010fa0378cp-20, 0x1.9b27169e51066p-21, -0x1.2cf2f079f069fp-23, -0x1.109cdd850f7a3p-21, 0x1.60e94dc9b6be1p-23, -0x1.502f4007e0e38p-27, -0x1.f28befe3d071dp-28, 0x1.aad1ea1a14a2fp-25, -0x1.e404a2820c0fdp-27, -0x1.5b4a90344b916p-28, 0x1.33eb20bb4e281p-29, -0x1.0140098b23c88p-28,
            0x1.5b6db6db6db9bp+0, 0x1.b6db6db6db6fep-1, 0x1.eb93e6122a0a4p-3, 0x1.08d669cf3f202p-2, -0x1.aa4497fd076afp-3, 0x1.e86bac1aee6b6p-7, 0x1.c4606a3f458f1p-5, -0x1.9b7654745ccacp-4, 0x1.012dc876cc4cep-5, 0x1.94a5764864a25p-5, -0x1.c690b97ddfffap-6, 0x1.aa949b4e58207p-9, 0x1.da6c83ba3d445p-9, -0x1.89a0114e8c0ebp-7, 0x1.11d924f871567p-8, 0x1.b287f6ac856f1p-9, -0x1.e3031104f9f14p-10, 0x1.2678db8df6ac3p-10, -0x1.76f1a5c2ab14cp-13, -0x1.ab58a8ef4ad89p-11, 0x1.4315a4b8659d1p-12, 0x1.fb86d71e436fap-16, -0x1.6618fa845ca6dp-15, 0x1.abc56858e285dp-14, -0x1.f07f0e8f12e4ap-16, -0x1.92e0df9b0ad1bp-16, 0x1.5c3aa4913d221p-17, -0x1.ec80746453a64p-18, 0x1.7cb4cdef85bddp-20, 0x1.02ea221b6a769p-18, -0x1.58f89f0d87144p-20, 0x1.102c380c1dadap-23, 0x1.f19f59e5c17adp-25, -0x1.a34670662945ep-22, 0x1.d40b40985ae94p-24, 0x1.7b9cfc0766bf3p-25, -0x1.53722116f2109p-26, 0x1.df3642512878cp-26, -0x1.bcc6af439ea92p-28, -0x1.f12b75535c9cfp-28, 0x1.4336914f6b5a6p-29, -0x1.6a79e1d8b53f0p-30, 0x1.d8e4ecbbe6c07p-33, 0x1.814e65af14a6cp-31, -0x1.abd2c967cf909p-33, 0x1.32e6d4ca5bdd2p-36, 0x1.9908c76682578p-38, -0x1.be07b708b61ffp-35,
            0x1.5b6db6db6db6ep+0, 0x1.db6db6db6db9ep+1, 0x1.0dba5cdc75457p-3, 0x1.73dac23d81853p-5, -0x1.4c412cd824e1bp-6, 0x1.b72d4ad82bca9p-6, -0x1.2eb699490c8c2p-8, -0x1.a695e963fa617p-8, 0x1.c5f919ca1f21ap-10, -0x1.39c97db9b19c1p-11, 0x1.7ffb14d7c2888p-16, 0x1.b22adabe5cbddp-12, -0x1.79582808b4f14p-14, -0x1.4510af572aa51p-16, 0x1.30b77b7ae07cbp-17, -0x1.2bc0ce61aa6ccp-16, 0x1.b08a68a773368p-19, 0x1.7b5c77c2ce146p-19, -0x1.8cd49a8d7c439p-21, 0x1.1499273b8fbe9p-21, -0x1.23fc6000fceeep-24, -0x1.7457eb213d09fp-23, 0x1.41f8e45b35aa5p-25, -0x1.7862eab9b3e22p-28, -0x1.42cf4cc7cd25ap-31, 0x1.03dfffb317306p-27, -0x1.89ad518fc979cp-30, -0x1.dda239051f7cbp-32, 0x1.365b02fa5ed1fp-33, -0x1.17d35f92da202p-32, 0x1.7326add2b6c56p-35, 0x1.412b2d1f5eab1p-35, -0x1.2683bf8190965p-37, 0x1.c4735d4422a6ep-38, -0x1.e947e7cb3f40ep-41, -0x1.f1c87ec1548cep-40, 0x1.8aed5964bcafap-42, -0x1.a89cb314b08ebp-44, 0x1.57835d4f362f9p-48, 0x1.2908916958edfp-44, -0x1.a9c01f0109fcdp-47, -0x1.8869d5253ca0ap-50, 0x1.662b4bb5da8d2p-51, -0x1.2355e28127c71p-49, 0x1.7c7a4acba9995p-52, 0x1.8aac9a7d39c89p-53, -0x1.6b71763fab2dbp-55, 0x1.d87c0f8ed2355p-55,
            0x1.6db6db6db6d83p+0, 0x1.3b6db6db6db78p+2, 0x1.aaf7749d8ea72p-4, 0x1.dda62336cf067p-6, -0x1.56ba8eb43895bp-7, 0x1.1ff80b0500f2cp-6, -0x1.5e47689d90b61p-9, -0x1.639ce37c707c9p-9, 0x1.3c667d8f7f6b9p-11, -0x1.71e916183cbe7p-12, 0x1.2c66a090dda3ep-15, 0x1.fd15f023e540ep-14, -0x1.7abc488b625aap-16, 0x1.35d6bc4be6eb5p-20, 0x1.7cb266e9ba796p-21, -0x1.072d15df31294p-18, 0x1.569413999f93ap-21, 0x1.246ed47f6fd2dp-22, -0x1.2ae9a55f71228p-24, 0x1.a036489be0225p-24, -0x1.d245da2032df8p-27, -0x1.02c3362dff9f4p-26, 0x1.95659c5018764p-29, -0x1.d3b483ec63a7cp-30, 0x1.840e129526966p-33, 0x1.2761561e304fap-31, -0x1.95d5bf4432c5ap-34, 0x1.7886d17ba389fp-37, 0x1.5ebc49017f7a2p-40, -0x1.08ac0d2a92d21p-36, 0x1.48a67190f3af5p-39, 0x1.75f2105c9573dp-41, -0x1.88c697ca6e476p-43, 0x1.8451964dfdbbfp-42, -0x1.b2ca722d7e70fp-45, -0x1.4eb86f1a4ebd7p-45, 0x1.fea0031347746p-48, -0x1.c8ec1f0baeb8ep-48, 0x1.b9ad8bb2c55eep-51, 0x1.69b6ed572b607p-50, -0x1.e46b7029a9bbap-53, 0x1.76cb2dbb93b32p-54, -0x1.e378b754b027bp-58, -0x1.35925aa2a0e45p-55, 0x1.7ba1bc6304fd5p-58, -0x1.ceb439c1d4d0cp-64, -0x1.12e2b367023edp-63, 0x1.bfe5d41415805p-61,
            0x1.9249249249247p+0, 0x1.2492492492486p-3, 0x1.fb19b8500895dp-4, 0x1.b018e998faee0p-2, -0x1.12face132cc43p-2, -0x1.dd9009f4dcc81p-3, 0x1.0f3967540fdbfp-2, -0x1.196320cb4c74cp-6, -0x1.99faaab80ff25p-4, 0x1.2f8fc98947e4dp-3, -0x1.61e97c07c8b1fp-5, -0x1.9c915eb12ad9dp-4, 0x1.f74aae75442edp-5, 0x1.a973ebfee77c7p-8, -0x1.1e407081dd5a5p-6, 0x1.be71a105c5170p-6, -0x1.19b42b69c0365p-7, -0x1.b64cc33e577d0p-7, 0x1.dc4884abee20ap-8, -0x1.47c65b34ef0d5p-10, -0x1.5dc44bc921923p-11, 0x1.a0971693338ffp-9, -0x1.2e4213cb955dcp-10, -0x1.7e9c303f0c205p-11, 0x1.c4bffa7aff6cbp-12, -0x1.63e95cfa0601dp-12, 0x1.29d0dc61a6b8fp-14, 0x1.a3f19a459d8b3p-13, -0x1.4c27dc54fc74cp-14, 0x1.2dad0bf18f236p-19, 0x1.01b3e4cae98a1p-17, -0x1.cd494fa3a56aep-16, 0x1.1abb541665e69p-17, 0x1.64ac456c70e1bp-18, -0x1.524268d7c0cd8p-19, 0x1.2cf840a528a67p-19, -0x1.05d3867d5bae3p-21, -0x1.0bc24a8f9191cp-20, 0x1.771b19b5be16fp-22, -0x1.20d7827157dfap-24, -0x1.00636dfc62d9bp-27, 0x1.dc2d93880d71bp-24, -0x1.14683e514e519p-25, -0x1.6be33d41c9308p-27, 0x1.71dee5bb2ab01p-28, -0x1.254ce4b76df76p-27, 0x1.192213c6381d0p-29, 0x1.2250e2cbab75dp-29,
            0x1.a4924924924c3p+0, 0x1.b6db6db6db6f0p-2, 0x1.3e763fe5d60bcp-3, 0x1.5c20ea2a404edp-2, -0x1.c19540862c80cp-3, -0x1.f37f6c04fc435p-4, 0x1.391a80ac6199dp-3, -0x1.751ba45b9dd30p-5, -0x1.1767e75e6e3bbp-5, 0x1.66673a03c0bd5p-4, -0x1.e62ff96e99cc7p-6, -0x1.5660edd668696p-5, 0x1.a9a1df93db4adp-6, -0x1.40ec5e1a77971p-9, -0x1.3302bc4495c35p-8, 0x1.7ad4df7435505p-7, -0x1.f31c0fdde0cedp-9, -0x1.1032cc9f9bb2bp-8, 0x1.2c2b26845cb63p-9, -0x1.970e446550bb2p-11, -0x1.0ba5137aece17p-14, 0x1.0315412191a2bp-10, -0x1.7174c3de2c344p-12, -0x1.57e7357a9a7efp-13, 0x1.b47ed7eec6282p-14, -0x1.b1179fa49e006p-14, 0x1.7f68bf713eedbp-16, 0x1.9414f49c0eef2p-15, -0x1.383b9232e4b73p-16, 0x1.27d1ccf360437p-19, 0x1.4c531aca63e55p-20, -0x1.a40e634f82648p-18, 0x1.f149e5af48f61p-20, 0x1.0e2c8efba3566p-20, -0x1.0260a36d57bb2p-21, 0x1.ff9f323d0a66dp-22, -0x1.b16166d51b8b3p-24, -0x1.8ff10703e4436p-23, 0x1.0eb037b05bbfep-24, -0x1.e5eb62a51e422p-27, -0x1.f74910874a95fp-31, 0x1.49f661bdc2a4bp-26, -0x1.6e3377c75269dp-28, -0x1.d3cf3e731e89fp-30, 0x1.d41e1503f1e00p-31, -0x1.71d138eb1f121p-30, 0x1.4f31e4e4d4566p-32, 0x1.6290d0027bdd3p-32,
            0x1.a492492492479p+0, 0x1.36db6db6db6a1p+1, 0x1.48086056aabf3p-3, 0x1.8f6dc96d0a8adp-4, -0x1.af1f044fcb851p-5, 0x1.eb504a49afef4p-6, -0x1.d0fa4d3d3009cp-11, -0x1.36db41a9a59f0p-6, 0x1.5dce07906d6b0p-8, 0x1.238f893fb8ce3p-9, -0x1.4509ca82c86c6p-10, 0x1.34ceaba6fc45bp-10, -0x1.70dae4f3c1247p-13, -0x1.e1b5e3d71d945p-12, 0x1.20d934b6b7fe1p-13, 0x1.5968616f865f3p-20, -0x1.c52852add39b9p-17, 0x1.2648bf7d732a7p-15, -0x1.fd2bbfb359abdp-18, -0x1.b6424dbb84250p-18, 0x1.2d8a714f013d6p-19, -0x1.61871199e13acp-20, 0x1.2bbe50118a753p-23, 0x1.5bae8e89f0069p-21, -0x1.60a7e772746d1p-23, -0x1.0c5ad25887d1dp-26, 0x1.0834e5bf7d26dp-26, -0x1.3f72289512867p-25, 0x1.011b780cc949cp-27, 0x1.a0321f4e28e89p-28, -0x1.fd0e1be910defp-30, 0x1.6401832ceff84p-30, -0x1.7bb2f0bced8fdp-33, -0x1.0e3b9e1e2373dp-31, 0x1.02e39b0fd69eep-33, -0x1.2c35cc88abd66p-37, -0x1.574a3e0b37c8ap-38, 0x1.c69cdccea2ef2p-36, -0x1.6f23299c4f4b1p-38, -0x1.592f3ff725ff0p-39, 0x1.bc08f4960de3fp-41, -0x1.0edb974dae6fcp-40, 0x1.64cea1352cb48p-43, 0x1.dcee7a1df000bp-43, -0x1.cab097a2d5bb1p-45, 0x1.89c0a4958a88cp-46, -0x1.1c6b80a17e726p-49, -0x1.96419c0ae7aa4p-47,
            0x1.a49249249249dp+0, 0x1.924924924926ap+2, 0x1.552619e49421ep-4, 0x1.5ca046b9bd8eap-6, -0x1.93c523992f465p-8, 0x1.75ee853aa6992p-7, -0x1.78c163b117e8bp-10, -0x1.549d7a6e9d205p-10, 0x1.f023addd83ab0p-13, -0x1.53b49421bb93ap-13, 0x1.fd549d21cc95ap-17, 0x1.48d62877e78aep-15, -0x1.95ad67f144e8ep-18, 0x1.d2442406ab51fp-21, 0x1.2e1972b107fb8p-24, -0x1.d6973710e93f2p-21, 0x1.01945b97978b3p-23, 0x1.3da34803ff78ep-25, -0x1.2c489f2d1baeep-27, 0x1.09aa871246eeep-26, -0x1.00215504775b3p-29, -0x1.c708cd0907b8dp-30, 0x1.32e63240d9faap-32, -0x1.c51617e5fa7ccp-33, 0x1.5f9b21c9789f5p-36, 0x1.7ee3149ff2f27p-35, -0x1.c1b9cd231eefep-38, 0x1.a16f744ef2325p-40, -0x1.c99aa0b915bb4p-46, -0x1.f5c7544786484p-41, 0x1.0b01cade8c08bp-43, 0x1.6ad154d3728ffp-46, -0x1.9b3bba67d9c28p-48, 0x1.0f3abbad82aa6p-46, -0x1.06493347160b8p-49, -0x1.445ec8760a9e0p-50, 0x1.b84fcdc695a25p-53, -0x1.df7fcb9610abcp-53, 0x1.982ff63536a36p-56, 0x1.1159882fa4315p-55, -0x1.3e9b4b371cc74p-58, 0x1.3957f00c54c56p-59, -0x1.8e7f4737447a9p-63, -0x1.6136afaed7bcdp-61, 0x1.771a246f77b24p-64, -0x1.1dc5ec832add3p-67, -0x1.fd733d2dcd7c8p-71, 0x1.7efd29765b026p-67,
            0x1.b6db6db6db6bfp+0, 0x1.0000000000005p+0, 0x1.6c5059873c12ap-3, 0x1.d9cbfe95a77cep-3, -0x1.2d7b7eed7ef6cp-3, -0x1.513a79bb453bap-6, 0x1.b9703118d3ab2p-5, -0x1.903471b6f909ap-5, 0x1.05311f17bcc7cp-8, 0x1.137e89ce38db6p-5, -0x1.bbd3c99f61f3fp-7, -0x1.a2e1551180162p-8, 0x1.50a60a9e3376cp-8, -0x1.de93544b52b44p-9, 0x1.7c06c34100471p-12, 0x1.484af41662560p-9, -0x1.e75bc36a635a8p-11, -0x1.2d10581c279a4p-12, 0x1.eeed2f07cd394p-13, -0x1.1be10546efbcdp-12, 0x1.d9fb8a55dbce3p-15, 0x1.e038f440a3033p-14, -0x1.6e5c4b1fc4467p-15, 0x1.e1d4c9a67803ap-19, 0x1.0be55e3493589p-18, -0x1.d59b504d3c7e1p-17, 0x1.fece8dba08c7bp-19, 0x1.6f16ea0055575p-19, -0x1.4158533c0cc74p-20, 0x1.c65a5e45a8d11p-21, -0x1.165ab425e3b80p-23, -0x1.ca7204d12a4ccp-22, 0x1.1c23c76164b8fp-23, 0x1.178284e881756p-29, -0x1.98c00a91bb8cdp-27, 0x1.38a1cf58c87d7p-25, -0x1.314b17e679f5dp-27, -0x1.a8ae157cc8e81p-28, 0x1.3fe5d93e4abd1p-29, -0x1.0252519e37058p-29, 0x1.5d03dfa27227fp-32, 0x1.99639cd826d72p-31, -0x1.d3be44f471ecdp-33, 0x1.f64cced5c6e0bp-36, 0x1.fbbdea09a4f0dp-38, -0x1.e6a5afa1676b8p-35, 0x1.cfd3ae74947efp-37, 0x1.79ea8473f7beep-38,
            0x1.edb6db6db6de9p+0, 0x1.2492492492477p-3, 0x1.d0c8b0fb833f7p-5, 0x1.5de5c64080124p-2, -0x1.f07c72bfeac8fp-4, -0x1.a57cf02b60d97p-3, 0x1.2e588dcf6a9f6p-3, 0x1.2913712e7248ep-4, -0x1.a22d5e84cc4f0p-4, 0x1.e9c21e051e5e6p-6, 0x1.b691809e597acp-6, -0x1.d97844703a9aap-5, 0x1.14a6079df3d8fp-6, 0x1.fd882a58ea171p-6, -0x1.2bc5831617100p-6, -0x1.9155bacb7fd38p-10, 0x1.5474b5a91d792p-8, -0x1.de180d04abb55p-8, 0x1.da3187918add5p-10, 0x1.e6de24a425bf2p-9, -0x1.daf0e570e18b0p-10, -0x1.d9617697dc74cp-16, 0x1.6090e63e06921p-12, -0x1.64a616a885da3p-11, 0x1.9683c7d4e5c01p-13, 0x1.ecad9dfee2803p-13, -0x1.d8ea00261364cp-14, 0x1.16e782d31fa36p-15, 0x1.24c6f99791549p-18, -0x1.6d26a6fb93ec2p-15, 0x1.d6555a9aea503p-17, 0x1.e01987118ca4bp-18, -0x1.0fc728bb5121fp-18, 0x1.f03ed398983c7p-19, -0x1.7d13d78ad1d79p-21, -0x1.cbe5f8ee1c37dp-20, 0x1.46cb93b694d72p-21, -0x1.26de215a95215p-25, -0x1.bfa66b88d4cb7p-25, 0x1.93d899be2efcdp-23, -0x1.b339607060dcbp-25, -0x1.222a3645b3067p-25, 0x1.f23937cb45296p-27, -0x1.8d713ba3a597ap-27, 0x1.13870373be83cp-29, 0x1.61852719d6823p-28, -0x1.b897db4a098c1p-30, 0x1.0fc0e94af5a22p-33,
            0x1.edb6db6db6da0p+0, 0x1.fffffffffffccp-2, 0x1.cf15155cd6e89p-4, 0x1.2bfed9d26c988p-2, -0x1.2518ac69a0b55p-3, -0x1.d5be537fa252ep-4, 0x1.b08f7e7ac3cc2p-4, -0x1.e758418cb4f8dp-13, -0x1.4b57bb74cf259p-5, 0x1.546606637887cp-5, -0x1.bf3acb1427e4bp-9, -0x1.e8e8689f4f98bp-6, 0x1.a4d62ce7493c9p-7, 0x1.d013dc4c8dd42p-8, -0x1.8530ed13277d4p-8, 0x1.aa867c4ab024ap-9, 0x1.1c2d714c7d4b3p-14, -0x1.85038e4f6d59fp-9, 0x1.1b5fcf82fe17bp-10, 0x1.3f4684f703119p-11, -0x1.ae9f3a31fae0dp-12, 0x1.23f8b0d497fdcp-12, -0x1.f67bc0d0e829dp-16, -0x1.8a87f5c603b36p-13, 0x1.1d7edd805fbf6p-14, 0x1.2f0b9f4f04e6cp-16, -0x1.021090e650f36p-16, 0x1.4d78049004f86p-16, -0x1.25e0bfd154a3bp-18, -0x1.f6d9a2f455478p-18, 0x1.8107317c3b0d4p-19, -0x1.f05a0d39ec548p-22, -0x1.8c6aac012a89dp-23, 0x1.f1841759dfda3p-21, -0x1.130e6543243b1p-22, -0x1.4ae2319a658a5p-23, 0x1.2ffa57d0b13a8p-24, -0x1.fac3c6addbaafp-25, 0x1.584223521400cp-27, 0x1.bbc0c30f06db9p-26, -0x1.17285788ed46bp-27, 0x1.e2bb8346c0558p-32, 0x1.42a7ced1c2c9dp-31, -0x1.36cee3e048d71p-29, 0x1.32eb0dae6c28ap-31, 0x1.7bfe7d01910aap-32, -0x1.2a069294961ebp-33, 0x1.07095fa0fb884p-33,
            0x1.edb6db6db6db5p+0, 0x1.9249249249264p+0, 0x1.4160c01f1dd41p-3, 0x1.4fd45ee004dbfp-3, -0x1.7045bedda2f41p-4, 0x1.621c2a25043cfp-9, 0x1.537e2320ad8f0p-6, -0x1.c925a3b7ee6f9p-6, 0x1.0a9cf2b5e04c4p-8, 0x1.952f3bfc0b0a8p-7, -0x1.2b5bad76f4546p-8, -0x1.3275abcb597d9p-10, 0x1.30855a2402579p-10, -0x1.3381c1cc4c156p-10, 0x1.52720dd75d322p-13, 0x1.18a6872b8ed51p-11, -0x1.7a83462c61fa8p-13, -0x1.0ea2ba925da09p-15, 0x1.1d38e70ddf8ccp-15, -0x1.8d4fa6bd58c86p-15, 0x1.2de5807158fd9p-17, 0x1.05a1ceb10646cp-16, -0x1.68569dbafa347p-18, 0x1.522ff67a28eb2p-21, 0x1.ca6655a372c51p-22, -0x1.9896ba9fcceebp-20, 0x1.81839ae8f1e91p-22, 0x1.248ad08006735p-22, -0x1.c8ba33cd6abd7p-24, 0x1.1dc25a42c8f65p-24, -0x1.e4f8993ba7cb6p-28, -0x1.1ebb5c5ae2e9cp-25, 0x1.35387c56001bap-27, 0x1.475e05fec8475p-30, -0x1.2190744844f81p-30, 0x1.230b940259069p-29, -0x1.d6df9eb58fe2dp-32, -0x1.d43518244bc27p-32, 0x1.27b1416d0f3e2p-33, -0x1.30a48ec7e44eep-34, 0x1.bd09f14468606p-38, 0x1.4252c88ba7940p-35, -0x1.396b4675fe99dp-37, -0x1.3f386dcb636f7p-40, 0x1.e67128e2cb8cdp-41, -0x1.1237dfb950bf4p-39, 0x1.abcd36a4fb2d3p-42, 0x1.69b6eba604d00p-42,
            0x1.2db6db6db6deap+1, 0x1.24924924924bep-2, 0x1.6780906438922p-5, 0x1.0ac24d6fbfec4p-2, -0x1.db538fd99d5aep-5, -0x1.ff2fb1cc2a59fp-4, 0x1.d4cdb8837d99dp-5, 0x1.9b9d285e19e87p-5, -0x1.55678999578ffp-5, -0x1.b449c48410d2ap-8, 0x1.405a8df38715cp-6, -0x1.730a89bb53471p-7, -0x1.91034630cc220p-9, 0x1.6bdf6b761d137p-7, -0x1.9093b2d48bfafp-9, -0x1.3157e82e359e1p-8, 0x1.4e7fd985de380p-9, 0x1.c1fdc31185002p-13, -0x1.7bd9b38fd455fp-11, 0x1.c0d278e704b71p-11, -0x1.1777b7a0a5306p-13, -0x1.d7f56d474a0abp-12, 0x1.7db0f353e7f88p-13, 0x1.979b1ea6180bap-15, -0x1.ade77dbbefbb8p-15, 0x1.b07c12c941cdbp-15, -0x1.03e6bfb6a8595p-17, -0x1.b2bc19314f90cp-16, 0x1.4096856e0a49dp-17, 0x1.e4a7db906afbfp-20, -0x1.11680838915ecp-19, 0x1.6ae60d24473c0p-19, -0x1.1ce94ad3bf9a4p-21, -0x1.0fe5c362b3a16p-20, 0x1.9275418bd1c40p-22, -0x1.618e18565d3bep-26, -0x1.6a7c8ea6e7137p-25, 0x1.e234ccb148380p-24, -0x1.d05fdf0ca98dap-26, -0x1.bade99ea9c4c4p-26, 0x1.646dba83bd5d6p-27, -0x1.3cbb8847e8e37p-28, 0x1.156dd89cf1728p-33, 0x1.cee0b1090c0c7p-29, -0x1.fc516e5a25d3fp-31, -0x1.606ccf43485dbp-32, 0x1.882c2f8153333p-33, -0x1.d56cd2689ec7dp-33,
            0x1.2db6db6db6de4p+1, 0x1.fffffffffffc4p-1, 0x1.ad5b8c28c6095p-4, 0x1.a4ea3d5f973a1p-3, -0x1.5474561d48d0ap-4, -0x1.9b00a7babc18dp-5, 0x1.4f481f100f2bdp-5, -0x1.0faccaa635fa8p-8, -0x1.7e4953e7c882ep-7, 0x1.9f6191c3a23f4p-7, -0x1.8a45b911d3639p-12, -0x1.ccfd1e5b6420ep-8, 0x1.2155ca2cc1699p-9, 0x1.c0cb328d0052fp-10, -0x1.10372cab2d782p-10, 0x1.148d3d3b05039p-12, 0x1.2dfc3df49dd79p-13, -0x1.840d838d18f0dp-12, 0x1.58e5e61abe55ap-14, 0x1.e726c79742520p-14, -0x1.9d3dcb5ce4f9ep-15, 0x1.a3c03c6554d0ap-19, 0x1.f0cfcce703605p-18, -0x1.f251e59a0aa73p-17, 0x1.9de756140a7a7p-19, 0x1.2aa48cc1dc40fp-18, -0x1.c952991022550p-20, 0x1.c94d9142da040p-23, 0x1.7a40b55ef80e0p-23, -0x1.0c4df14a6b471p-21, 0x1.e3eefef85bd53p-24, 0x1.e40c952c86f61p-24, -0x1.78d783e2952b1p-25, 0x1.16261c3311317p-26, 0x1.dbd1a5f73b6b7p-31, -0x1.cb8991902fa9cp-27, 0x1.cc1bea6630429p-29, 0x1.d3ea08637dda1p-30, -0x1.aee826a53db8cp-31, 0x1.6f5d5120d4b99p-31, -0x1.8125c40baac7cp-34, -0x1.2099b7dda6bdcp-32, 0x1.4011ef776ce2cp-34, 0x1.0bb6dc6841eb9p-38, -0x1.0ca02bcdd94eap-37, 0x1.2ed0debb247f2p-36, -0x1.e0a1f118e649bp-39, -0x1.d61b323ac9ae9p-39,
            0x1.3fffffffffffdp+1, 0x1.7fffffffffff0p+1, 0x1.ce5924b880e0ep-4, 0x1.699b80ce6e66bp-4, -0x1.1c316587dc011p-5, 0x1.3db77f3b71755p-7, 0x1.812d3cf3957dap-9, -0x1.14d60229bed67p-7, 0x1.5d7c32dfcb193p-10, 0x1.c7f32d213d6f9p-10, -0x1.143771b2ee6e8p-11, 0x1.a5b87db7bd480p-15, 0x1.c41ec0746f2b9p-15, -0x1.dd42b5c7aa1cfp-14, 0x1.1bbbfa3e8d8e7p-16, 0x1.957b33e0da855p-16, -0x1.c5dd87de80c6ep-18, 0x1.93883d2c7b595p-21, 0x1.2e1a034573ef8p-21, -0x1.826ecb05ce502p-20, 0x1.ef2cd974f076ep-23, 0x1.1496f8f1dca0fp-22, -0x1.3ab8ef9f98754p-24, 0x1.64c7e4f0832dap-26, 0x1.749a8ac2e08ebp-29, -0x1.1821774513bc4p-26, 0x1.92bebc6dab480p-29, 0x1.0fc26c442299dp-29, -0x1.5b314cf1a097ep-31, 0x1.b2f7f3f2cbc0fp-32, -0x1.9deca1f8c017fp-36, -0x1.50ab4661cd542p-33, 0x1.0efe92cc2e181p-35, 0x1.0ed781b39d89dp-37, -0x1.0b647cbb5d08fp-38, 0x1.6c6ab0126894ap-38, -0x1.7e540d7def77fp-41, -0x1.3223e46ae651dp-40, 0x1.17950a048590fp-42, -0x1.a67c3c770bd4dp-45, -0x1.385b385e923b0p-47, 0x1.b41399b1e1ecbp-45, -0x1.254833215e268p-47, -0x1.5ac0479895d89p-48, 0x1.9054cccd09c8ap-50, -0x1.4865ada20dad1p-50, 0x1.05fd801a5141bp-53, 0x1.6e2ee0c0218a1p-52,
            0x1.6492492492455p+1, 0x1.fffffffffffcep+0, 0x1.a3b4f40da2a23p-4, 0x1.0ba8653fbffe1p-3, -0x1.8a357009a683ep-5, -0x1.3babd3892db9dp-7, 0x1.939d4fc7859f4p-7, -0x1.f2f855f398fb2p-8, -0x1.efe412bf7656ep-11, 0x1.0f3f8f9800ed6p-8, -0x1.703979a853188p-11, -0x1.095ad875aa432p-10, 0x1.703118fe2a33fp-12, 0x1.1f286d9ac4a5fp-15, -0x1.21eb0c41ae52ap-14, 0x1.2c0ecff6adf75p-14, -0x1.293ea56c9ed30p-18, -0x1.c8b4d80d20a4bp-16, 0x1.c544fd45d4578p-18, 0x1.b949fe005290cp-19, -0x1.d1e83b48025e0p-20, 0x1.0facbe5ccd469p-20, 0x1.6864fb4978a34p-26, -0x1.1b0761440051ap-21, 0x1.e660149722c19p-24, 0x1.43ba74b5bf08fp-24, -0x1.09d923657d3dfp-25, 0x1.03c6e1c989009p-26, 0x1.0fc97c732d8e2p-31, -0x1.25e72b3c8883cp-27, 0x1.e9729a77d05a8p-30, 0x1.2f16b7b4977edp-30, -0x1.dc5c7d3e6a723p-32, 0x1.199af90603838p-32, -0x1.0e4b4727b56cbp-37, -0x1.0506afa74949dp-33, 0x1.c1a0f0892bc02p-36, 0x1.742de232b9ab3p-37, -0x1.4b346c5192ad4p-38, 0x1.2a06b96011603p-38, -0x1.d746b7b5365b9p-42, -0x1.8011c7493afb7p-40, 0x1.61f4c85efc035p-42, 0x1.6b02b45d68ed9p-45, -0x1.430a9ce9d27ccp-45, 0x1.0b69cf3177ee1p-44, -0x1.397d34c19554ep-47, -0x1.b6d648040857fp-47,
            0x1.6db6db6db6da8p+1, 0x1.6db6db6db6de0p-2, 0x1.f9a8af0cf6a76p-6, 0x1.a91430d6b52adp-3, -0x1.cc751a5a9aaf6p-6, -0x1.467fd895827c2p-4, 0x1.57818b697695cp-6, 0x1.ef3d64ce3c8efp-6, -0x1.cad29b616aef4p-7, -0x1.34df3f4b11e6ap-7, 0x1.00cc998b281dep-7, 0x1.bf048378b8947p-11, -0x1.a7dbaaaa97fa3p-9, 0x1.b777956b4a2f2p-10, 0x1.5dbc933729525p-11, -0x1.881fdd1ad3823p-10, 0x1.db89c55bf61b7p-13, 0x1.616448805ef34p-11, -0x1.19990325802b4p-12, -0x1.0be6868f0bf2dp-13, 0x1.cd4b6a27ca082p-14, -0x1.8d545cd79ccd0p-15, -0x1.a4129b36b1554p-17, 0x1.777774470bc78p-15, -0x1.633b0701092d4p-17, -0x1.cff22c53aedd2p-17, 0x1.a71b7e46351fcp-18, -0x1.07abf3992befbp-22, -0x1.3657104e8e7cap-20, 0x1.fafea619a08aep-20, -0x1.64f4e41afd183p-22, -0x1.6b5f481eb2adbp-21, 0x1.0ae42cd70b8b0p-22, 0x1.774ac45854496p-26, -0x1.8d36830e06916p-25, 0x1.1a346e94fb606p-24, -0x1.77d7388a4ae7ap-27, -0x1.8635330dc29aep-26, 0x1.08e27cf827b23p-27, 0x1.7a2171397904fp-32, -0x1.3d43ec1f1af4ep-30, 0x1.128c5f02e8468p-29, -0x1.96998b795faf6p-32, -0x1.3603a50975cedp-31, 0x1.a4d371306957ap-33, -0x1.aa26fecbf74a3p-36, -0x1.28dc9cf6849e4p-36, 0x1.c39f1d0ab07dbp-35,
            0x1.8000000000034p+1, 0x1.1249249249261p+0, 0x1.187e15c2a5a27p-4, 0x1.5daa74024cfc7p-3, -0x1.70c647552fc80p-5, -0x1.5e689cf434b18p-5, 0x1.53cd2932606f6p-6, 0x1.7d6ce32c58b30p-8, -0x1.df2ac957acfe7p-8, 0x1.0e44bcbd22057p-9, 0x1.be2a38b365e45p-10, -0x1.116f4b20f6f45p-9, -0x1.8d336e23e73edp-17, 0x1.f1a14bdd2996ep-11, -0x1.c60c49464cafcp-13, -0x1.05aaf34c4de53p-12, 0x1.d824eb3918c44p-14, 0x1.2d8f8526d7f4fp-17, -0x1.d3f2f66e0f7a8p-16, 0x1.abb7918579c2fp-16, -0x1.7856db08e1259p-22, -0x1.9e7851e3a6694p-17, 0x1.a5af182e8c60ep-19, 0x1.407612062e0f8p-19, -0x1.3e0eb99c0548ep-20, 0x1.63ea72157780fp-22, 0x1.12da303892a9fp-23, -0x1.7a0e3a35d0d22p-22, 0x1.17eeffd6c3aa3p-24, 0x1.8450209d836cep-24, -0x1.16215559d34b4p-25, 0x1.b0f48039e4e6bp-30, 0x1.3449a5e4c12cap-28, -0x1.14eaaf4941a20p-27, 0x1.61648e7889a49p-30, 0x1.31bd04c3d98dfp-29, -0x1.8173fc2623e80p-31, 0x1.114dabdfd247fp-40, 0x1.90211d6f7c455p-34, -0x1.68b28f5b1310fp-33, 0x1.d38014aef32b4p-36, 0x1.6c76d11fc6525p-35, -0x1.bc679e72b645ap-37, 0x1.3ca4ac3f6e23ap-40, 0x1.4fb61f45ffe1bp-40, -0x1.a3380de4588bcp-39, 0x1.2bf5a4c74a6b6p-41, 0x1.4e5c2596caf1ep-41,
            0x1.800000000003fp+1, 0x1.0db6db6db6db8p+2, 0x1.6dd4ed9682dd5p-4, 0x1.f63d03625a0dbp-5, -0x1.39a89f593279bp-6, 0x1.fafe14b875818p-8, 0x1.79a8d557d5416p-11, -0x1.f6fd6d024e8c1p-9, 0x1.1b6ad3b179104p-11, 0x1.02ea89ca0006fp-11, -0x1.11878531ee718p-13, 0x1.467a8b1f039a1p-15, 0x1.8b5c65db84421p-18, -0x1.98efce4654de3p-16, 0x1.c0e6c5bd4d3f3p-19, 0x1.a37f929e696ccp-19, -0x1.a654f669b39cdp-21, 0x1.150748d759f92p-22, 0x1.b906e4ffae318p-26, -0x1.39b5726cee97bp-23, 0x1.66e7a95e3837ap-26, 0x1.1319233b62baap-26, -0x1.1eb82eabf1567p-28, 0x1.12f4d405dacbbp-29, 0x1.6c14370366f2bp-37, -0x1.baec4b9ebed6ap-31, 0x1.1163e4437e6abp-33, 0x1.114afd5681923p-34, -0x1.48afe83db2153p-36, 0x1.01fabd1afcc21p-36, -0x1.04410250f3e79p-40, -0x1.136343110d47ap-38, 0x1.74ab0b09f928ep-41, 0x1.24046cd15a9f8p-43, -0x1.21d5e0b1bb17fp-44, 0x1.9d8f1865f0cd3p-44, -0x1.5b4737601f343p-47, -0x1.1be354ee06e66p-46, 0x1.ae1c32b963eadp-49, -0x1.f784542693586p-52, -0x1.185becfe153c0p-53, 0x1.108dd72245648p-51, -0x1.213d8f4541b8ap-54, -0x1.b0d735d203ebfp-55, 0x1.896f2b2e593b6p-57, -0x1.dbc74227d31f0p-58, 0x1.7d58d358c9335p-62, 0x1.1e2d3a2b647d5p-59,
            0x1.c000000000031p+1, 0x1.b6db6db6db6c8p-2, 0x1.6e99573642bb7p-6, 0x1.531cb6b4a4901p-3, -0x1.e1b469011eb86p-7, -0x1.9892809d2ee62p-5, 0x1.f8b7c93874909p-8, 0x1.f4a9812a30537p-7, -0x1.f29e794e292c5p-9, -0x1.2d9384cba73c2p-8, 0x1.dd3ca44d85654p-10, 0x1.42143e9c99e5cp-10, -0x1.ac0bb2a9f796bp-11, -0x1.919b26af16133p-13, 0x1.4dc573cf94643p-12, -0x1.237aad93170fcp-14, -0x1.8862e8db81c60p-14, 0x1.657d12caefe6dp-14, 0x1.7ef3f431db8dbp-17, -0x1.8bc85833b35d0p-15, 0x1.f1e5f6c8489d8p-18, 0x1.17a8895084882p-16, -0x1.95d3be6e97bcbp-18, -0x1.80efb93fcaa5cp-19, 0x1.317551acaae2dp-19, -0x1.8a4823c8935a5p-21, -0x1.8f979ed260f89p-22, 0x1.8f5923b44ab68p-21, -0x1.b07d373eb6450p-24, -0x1.1a672fe9ea19fp-22, 0x1.776c4de8bc7cdp-24, 0x1.0dc6361db0b76p-25, -0x1.b74c2983d963ap-26, 0x1.08dc0e652844cp-26, 0x1.2ab87687306b9p-30, -0x1.3f925a1016981p-27, 0x1.1f80f60d6b406p-29, 0x1.0d8b00934fdaep-29, -0x1.ce44d5ee38bf4p-31, 0x1.784d90a2cfcffp-33, 0x1.abfec3e2fc149p-34, -0x1.fd3028acb5925p-33, 0x1.60be845dea4a9p-35, 0x1.06d4146834cd1p-34, -0x1.634c3e03324e0p-36, 0x1.17061680e0423p-41, 0x1.8a079cb4520a8p-39, -0x1.5205b911cf8b6p-38,
            0x1.c00000000001cp+1, 0x1.56db6db6db6a3p+2, 0x1.2ee6be1f13165p-4, 0x1.826f0368f722cp-5, -0x1.9108725e3bfcep-7, 0x1.7f3aabd24e9a1p-8, 0x1.ba565122f1f73p-13, -0x1.0ffa7b23259a8p-9, 0x1.0c7f772cfc823p-12, 0x1.8e87bec776b30p-13, -0x1.75ad986296eedp-15, 0x1.4024af7d7ed77p-16, 0x1.bdde99bb6c0efp-21, -0x1.e81450f93b714p-18, 0x1.df0e6261d736ap-21, 0x1.58599156027b6p-21, -0x1.3d8680b775acfp-23, 0x1.30162b403e0edp-24, 0x1.3437e3b5cbaa3p-30, -0x1.a10b513eb741ep-26, 0x1.a8345ca36f59dp-29, 0x1.f22b74b6b0452p-30, -0x1.e71b5b0717dc8p-32, 0x1.33406d43f6daep-32, -0x1.f1b39c7438f5dp-38, -0x1.4d7f50267c01bp-34, 0x1.68f148e59a35bp-37, 0x1.17aabc8c4a971p-38, -0x1.46928de8cb59bp-40, 0x1.2d772de71d961p-40, -0x1.2818dcc3bc3f3p-44, -0x1.e4ea915d80b5cp-43, 0x1.1c4214bb1d8c3p-45, 0x1.4f5c0ffc9399fp-48, -0x1.6446ad3805ed9p-49, 0x1.0e224643f6517p-48, -0x1.7ef041a029448p-52, -0x1.334f7573fd79cp-51, 0x1.8f0f0852eb84fp-54, -0x1.a8ef910d4f4e5p-57, -0x1.f79a06e7a8ae6p-59, 0x1.aa9789512b160p-57, -0x1.753417234a9a3p-60, -0x1.3aef816413228p-60, 0x1.dd26a78d650b9p-63, -0x1.da24bcf44806ap-64, 0x1.09258ecab3626p-69, 0x1.206392faed2b0p-65,
            0x1.d24924924922dp+1, 0x1.6db6db6db6dd8p+0, 0x1.d5ce949c3e962p-5, 0x1.16a68dbc44f23p-3, -0x1.dc9a16f0d5896p-6, -0x1.b6989e74612dep-6, 0x1.54068ff895f92p-7, 0x1.8b2fa0eec9990p-9, -0x1.80dc93b7baeb1p-9, 0x1.506cd0510b829p-11, 0x1.441724afcecdbp-11, -0x1.2f8aa988d2132p-11, -0x1.deb76c94e6d51p-15, 0x1.ea27f1d6c006dp-13, -0x1.ab1b2349dae63p-16, -0x1.0af454ff52ca4p-14, 0x1.1a3f4af16d9dap-16, 0x1.4abfb38ff62f6p-17, -0x1.6cd4c8e35ffa2p-18, 0x1.f1603383ba66bp-21, 0x1.00c2cd497d56fp-20, -0x1.4609fb948530bp-20, 0x1.748d6fc4dacd9p-25, 0x1.c9c7862f7720bp-22, -0x1.8ab428fb469a5p-24, -0x1.38b9783a602f3p-24, 0x1.0c464d80b436ep-25, -0x1.5573a836bd497p-28, -0x1.324791bb2ad22p-28, 0x1.e9746e7a7f3f0p-28, -0x1.7a9043cab9682p-31, -0x1.1edce2df8b3cbp-29, 0x1.24027fabad7c3p-31, 0x1.ad6c41d7d6b6bp-33, -0x1.09e916515d81ep-33, 0x1.533686db36234p-34, 0x1.433ca8740e80fp-39, -0x1.3a126e0f56b9fp-35, 0x1.e138913f79e01p-38, 0x1.844bc4bb21114p-38, -0x1.1e4d0670684a4p-39, 0x1.59c703a32d900p-41, 0x1.4667b5e6c7fd9p-43, -0x1.10ff7f9a0e89ap-41, 0x1.5ce8fc9698cabp-44, 0x1.957112071a164p-44, -0x1.e8f58707e7382p-46, 0x1.2f61d17a3494ap-48,
            0x1.ffffffffffff7p+1, 0x1.5b6db6db6db3ap+1, 0x1.16747a225544dp-4, 0x1.88b1df9ad9b09p-4, -0x1.7f78df2831e5ap-6, -0x1.f548cbc4a0dc4p-8, 0x1.3fa54bda5e988p-8, -0x1.ca454f087a58ap-10, -0x1.28df5171a5059p-11, 0x1.ea60141d81049p-11, -0x1.354281501453fp-15, -0x1.de53aaaa55d32p-13, 0x1.3d7e595dc4b28p-15, 0x1.02a84feae45a2p-15, -0x1.58dbde373d767p-17, 0x1.68bacc20f488dp-21, 0x1.861266979f44bp-20, -0x1.af778ea250679p-20, 0x1.b75667e4b4b2dp-26, 0x1.e6bc996f84cc0p-22, -0x1.3634322c8bccep-24, -0x1.06c387578850cp-24, 0x1.4e364ff33168dp-26, -0x1.5905b80fc3373p-29, -0x1.35f489bd69b03p-29, 0x1.cc3eeff229e0bp-29, -0x1.ec5df2dc3e9acp-33, -0x1.b46f9e30b4653p-31, 0x1.50f842bc792adp-33, 0x1.25fd6af9e6a21p-34, -0x1.04df32f6fa688p-35, 0x1.15dd781c1a2afp-36, 0x1.363bdff6ec03ap-40, -0x1.d0abd48d53c04p-38, 0x1.01becdb2efb39p-40, 0x1.0c5ab177d0ce0p-40, -0x1.20636443a457dp-42, 0x1.4288e97361033p-45, 0x1.9c7335670e4cfp-46, -0x1.91b4157b1b65bp-45, 0x1.49fe7c83a0175p-48, 0x1.311035a73ca5ep-47, -0x1.06be62b5b7647p-49, -0x1.883f3697d7877p-53, 0x1.ee64562bc0d12p-53, -0x1.37a4584d18057p-52, 0x1.84c578719e6a1p-56, 0x1.18a120ab58b12p-54,
            0x1.124924924922cp+2, 0x1.0000000000008p-1, 0x1.0f08ebb8f2e57p-6, 0x1.114b01291e7edp-3, -0x1.119fb662ee243p-7, -0x1.05e16ad90ec43p-5, 0x1.ab5f38eb52b33p-9, 0x1.fb0272e3617a2p-8, -0x1.32c9101a48109p-10, -0x1.ec865f88622edp-10, 0x1.ab938ff4bca48p-12, 0x1.d859d06baa47fp-12, -0x1.275299dddfac5p-13, -0x1.ace10402b21a7p-14, 0x1.9465066651c45p-15, 0x1.454bc972ed121p-16, -0x1.0cdc2fef05240p-16, -0x1.7b7329c6e18b4p-20, 0x1.4b9c4be5848fep-18, -0x1.64e9192adcecbp-20, -0x1.5afce735a5d69p-20, 0x1.14dac2918ca58p-20, 0x1.de8596e5913cap-23, -0x1.06722b2e1878ep-21, 0x1.2334796a5f7eep-26, 0x1.73d05ddc4a945p-23, -0x1.2b140f265de04p-25, -0x1.769a673a57c6cp-25, 0x1.23f4ff7bb6cb9p-26, 0x1.49821dd9641adp-28, -0x1.5da81a4f06accp-28, 0x1.1dec741a588c3p-29, 0x1.acfe504a7f7e3p-31, -0x1.99440ecfdf809p-30, 0x1.12b8da5919c57p-33, 0x1.0cb9a7a85211fp-31, -0x1.152a45b435ad1p-33, -0x1.59df00ba44394p-34, 0x1.6f679d8cbcd7bp-35, -0x1.5af62e9c81be9p-37, -0x1.ace2ad2e570c5p-38, 0x1.726f80fbb9efap-37, -0x1.2449ec242c41bp-40, -0x1.c55bf3aa77283p-39, 0x1.df04d0a4b83c9p-41, 0x1.971587cd98fd7p-42, -0x1.f0f1081633805p-43, 0x1.fe41f52079065p-44,
            0x1.36db6db6db6fcp+2, 0x1.0492492492489p+2, 0x1.dc32ab83bc239p-5, 0x1.1503fa71d7982p-4, -0x1.ca1b72748e2ccp-7, -0x1.fe43e0c43db3ap-10, 0x1.e4e5498839e48p-10, -0x1.495b7f9f35da6p-10, -0x1.4deb550f5e063p-14, 0x1.6a16a1d888254p-12, -0x1.df231609c6760p-16, -0x1.959d0673d4461p-15, 0x1.200cc6c3b35b5p-17, 0x1.385397d3d927fp-19, -0x1.47ca82ed8aef8p-20, 0x1.7799eebdbe9bcp-21, 0x1.de96ca306e08dp-25, -0x1.ce49018c4447fp-23, 0x1.3e0c12a3c0e2dp-26, 0x1.ffb07d820bb37p-26, -0x1.74270e2602bc1p-28, -0x1.23760c19309d8p-30, 0x1.7dfd6dedaa8b3p-31, -0x1.1c3ca0f32b153p-31, -0x1.9c0903dc02729p-37, 0x1.28e9e6b93c647p-33, -0x1.02afbfa6d7072p-36, -0x1.0b246098a6a6cp-36, 0x1.cd5902ee6cb37p-39, -0x1.e337b75798e85p-43, -0x1.50581dac50d25p-42, 0x1.cf68ee471a0dbp-42, -0x1.6b7478cb2adaap-46, -0x1.553f3001359dcp-44, 0x1.8d5a583553e8cp-47, 0x1.5b949f276204bp-48, -0x1.d25e044e77795p-50, 0x1.f893a910c61cep-51, 0x1.93c7d7559bc1bp-55, -0x1.3ea412c611822p-52, 0x1.0ffe715dde5b9p-55, 0x1.187367f89a5c0p-55, -0x1.d27262962066dp-58, 0x1.fd0446ba94ebbp-61, 0x1.09ad3a45b2ffcp-61, -0x1.f1f0dd4ed12bfp-61, 0x1.2cf4b1e6f6c18p-64, 0x1.2d0b242ff676ap-63,
            0x1.3ffffffffffdep+2, 0x1.a49249249245ep+0, 0x1.1ff6f059a5bf6p-5, 0x1.a5b545715209ep-4, -0x1.b1af913a9e11bp-7, -0x1.114c566536d3fp-6, 0x1.da5d0810e3911p-9, 0x1.1f3442de5e401p-9, -0x1.b8a2bf958ebd2p-11, -0x1.34c3355d1ffafp-13, 0x1.65d2e1c3f7e86p-13, -0x1.37b1dde8a22e5p-15, -0x1.eeb2ded255073p-16, 0x1.61b04c5145837p-16, 0x1.ffe70cc6dcf23p-19, -0x1.bb1cd121adcd7p-18, -0x1.0885be2c8d354p-23, 0x1.b2a5d472e395cp-20, -0x1.20d426c5b2f9cp-23, -0x1.5a778f3c4b608p-22, 0x1.fae9189d05f4ep-25, 0x1.9aa8f1ae44b68p-25, -0x1.19d5371d7e8d3p-26, -0x1.6de775581404dp-29, 0x1.d1bfae7f09540p-29, -0x1.79128a1bf0b48p-30, -0x1.ffc22c022469fp-32, 0x1.6a48caa132c7ap-31, 0x1.c82d27589b9f1p-39, -0x1.87f0d97280829p-33, 0x1.891d60f8f0cf8p-36, 0x1.1d0fa0f9cdb09p-35, -0x1.2197536de75f4p-37, -0x1.80f9f3926e145p-39, 0x1.f79c5b29c847fp-40, -0x1.6e88a70986b91p-41, -0x1.d8da592937475p-43, 0x1.94336794fd8c8p-42, -0x1.2b23fe3cc24b9p-46, -0x1.8fc98e9d39c66p-44, 0x1.1716f8d4010fbp-46, 0x1.a9423e54e9faep-47, -0x1.2acb9c0d6783dp-48, 0x1.effc7d22ff349p-52, 0x1.4682ad57bf9d4p-51, -0x1.881a5d255af42p-51, 0x1.4a97cb1eb2342p-56, 0x1.abad4f9ebeaafp-53,
            0x1.56db6db6db6dap+2, 0x1.b6db6db6db6f1p-2, 0x1.21d18783b4413p-7, 0x1.b43f860dbfab7p-4, -0x1.c720e525a3d14p-9, -0x1.4d31459fd15fap-6, 0x1.1114a22dc1fbep-10, 0x1.01655b0cc13e9p-8, -0x1.2921e3723f3a0p-12, -0x1.921bf5bb38246p-11, 0x1.35961762354d1p-14, 0x1.3d521e163c225p-13, -0x1.3cc44e9940d31p-16, -0x1.f9130c46d9239p-16, 0x1.42e11b67bcc55p-18, 0x1.93edd5f6becbap-18, -0x1.4ae479a3ba064p-20, -0x1.4244b496ec9c9p-20, 0x1.5701877b30a3cp-22, 0x1.f910b5276ba37p-23, -0x1.68a6fb4a2d702p-24, -0x1.76b96c560cf50p-25, 0x1.7fcb86c983f2fp-26, 0x1.da1d2160f52dap-28, -0x1.9930ed068f5e8p-28, -0x1.22f5d1f00e21fp-31, 0x1.ab8d505a0d605p-30, -0x1.1515a8916f316p-32, -0x1.a48ab06e45dc0p-32, 0x1.9a436cab29736p-33, 0x1.6858e9e70ff2bp-34, -0x1.6ea1597dbfc68p-34, -0x1.b06b2044a38d8p-37, 0x1.07c4918e67bf9p-35, -0x1.c84f7ec061a31p-42, -0x1.3cb9b86096809p-37, 0x1.5f1ffbdead253p-40, 0x1.30d6dc21a9aeap-39, -0x1.4e9735cb4aeb5p-41, -0x1.885c9415d4690p-42, 0x1.a7a16d88419cbp-43, -0x1.50bdc98774bc5p-48, -0x1.763d1c9401b1ap-45, 0x1.07d10a8ed55a9p-45, 0x1.4295a42f68c16p-48, -0x1.cd44020f54c34p-47, 0x1.3c830c00f68acp-50, 0x1.e9b32707eb9cdp-49,
            0x1.8000000000036p+2, 0x1.5249249249252p+1, 0x1.243076e59a7d3p-5, 0x1.43cdf3367d3f8p-4, -0x1.4cc867c3b560bp-7, -0x1.20383b077222bp-7, 0x1.08dd92d53534dp-9, 0x1.26123d96b6451p-11, -0x1.541418ef1f96fp-12, 0x1.f62fa8fdf0d1cp-15, 0x1.5dca643872a25p-15, -0x1.08034dcfa4693p-15, -0x1.e5876a79db31ep-19, 0x1.f7d2b2c33f9a5p-18, -0x1.b36f5ce79a893p-25, -0x1.66fda8b4d1ccdp-20, 0x1.db2807885ce35p-24, 0x1.8f1edc16a6115p-23, -0x1.05cab0910b771p-25, -0x1.2d5cfe3d58883p-26, 0x1.8b093c9a3fa11p-28, -0x1.d40aa5e0b24dep-37, -0x1.bd47049678174p-31, 0x1.0b0ff461b7d90p-31, 0x1.3be73388a369fp-34, -0x1.302fb2ce89387p-33, 0x1.3a720fbb2496fp-39, 0x1.c49124499e8c0p-36, -0x1.75630c585640dp-39, -0x1.dad1b06676d6dp-39, 0x1.7f6d40bb7981dp-41, 0x1.ebdb3b21a635fp-43, -0x1.01160b3f2702ep-43, 0x1.3c0d77b43c65cp-45, 0x1.b3429623a1240p-47, -0x1.1e7e50f6d0934p-46, -0x1.78a0f91149b1fp-55, 0x1.e1973c3709bcep-49, -0x1.80d983de7c5cfp-52, -0x1.02afb91503b1ap-51, 0x1.a059e8af84c90p-54, 0x1.c466e0d3818e1p-56, -0x1.07eb8a17ff8e1p-56, 0x1.d0d9d1327daa3p-58, 0x1.573ca65f8c63ap-60, -0x1.5883ac5de70bdp-59, 0x1.f1e2543545ecep-64, 0x1.f1531951949fep-62,
            0x1.a49249249246bp+2, 0x1.b6db6db6db6b3p-1, 0x1.7591fbe3ea040p-7, 0x1.5d90b4c7c4375p-4, -0x1.cef86ada94f3bp-9, -0x1.a3f7b6a27032ap-7, 0x1.b0d3baecf1570p-11, 0x1.f3c906981ff6dp-10, -0x1.69b22158fb150p-13, -0x1.257de0386417fp-12, 0x1.1cd7b98d26863p-15, 0x1.521ffa05638d4p-15, -0x1.b0b31fc6001dep-18, -0x1.7a2dd1582e205p-18, 0x1.40bf475be8297p-20, 0x1.92d53511b80eap-21, -0x1.d307bf3b9b545p-23, -0x1.883f0c796cbb2p-24, 0x1.4ef9162b43a22p-25, 0x1.372874c828873p-27, -0x1.d9b1e00550738p-28, -0x1.88f8364f12541p-32, 0x1.499cd64ccd04fp-30, -0x1.5f51da3999ac9p-33, -0x1.c14aade6135b6p-33, 0x1.2f6c0d501ceb0p-34, 0x1.2915622417950p-35, -0x1.5b600778f0fb6p-36, -0x1.761e5d7e0e4f9p-38, 0x1.587a4f95b793ep-38, 0x1.aef6663f4cda4p-41, -0x1.3c4ea2d5932d0p-40, -0x1.96c0d8a829f6ap-44, 0x1.13308f7c87da0p-42, 0x1.51840127d8b1fp-48, -0x1.c7ffd09370259p-45, 0x1.11b4267ac5085p-49, 0x1.651f8908c6ba3p-47, -0x1.159b7f0772c89p-50, -0x1.026d96ee73558p-49, 0x1.5f18478bb2207p-52, 0x1.47589db60e153p-52, -0x1.6e7919d4682cfp-54, -0x1.33b76db1c4bd7p-55, 0x1.4c30234131854p-56, 0x1.814f1a40b3884p-62, -0x1.02df2c77f8831p-58, 0x1.a77900467613fp-60,
            };
            // clang-format on
        } // namespace taylor
    } // namespace faddeeva_impl
} // namespace eos
