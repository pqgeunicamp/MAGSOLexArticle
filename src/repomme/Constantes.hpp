#ifndef CONSTANTES_HPP
#define CONSTANTES_HPP

#include <vector>				// cria e manipula vetores
#include <Eigen/Dense>

namespace Constantes
{
    //Informações das moléculas
    struct Molecula
    {
        //representação
        std::vector < std::vector <int> > SOLex;  
        std::vector < std::vector <int> > MAG;              // antigo Ligas
        std::vector <int> MAGnode;                          // antigo nElementos
        std::vector < std::vector <int> > LigacoesPorAtomo;
        
        int UltimoC;

        double Tm;      //Temperatura de fusão [K] calculado por MG
        double Tb;      //Temperatura de ebulição [K] utilizada da molécula i
        double Tb_IFP;  //Temperatura de ebulição [K] calculado por IFP
        double Tb_MG;   //Temperatura de ebulição [K] calculado por MG
        double Tc;      //Temperatura crítica [K] calculado por MG
        double Pc;      //Pressão crítica [bar] calculado por MG
        double Vc;      //Volume crítico [cm³/mol] calculado por MG
        double Gf;      //Energia livre de Gibbs padrão de formação @298K [kJ/mol] calculado por MG
        double Hf;      //Entalpia padrão de formação @298K [kJ/mol] calculado por MG
        double Hv;      //Entalpia padrão de vaporização @298K [kJ/mol] calculado por MG
        double Hfus;    //Entalpia padrão de fusão [kJ/mol] calculado por MG
        double LogKow;
        double Fp;
        double TAit1;
        double TAit2;
        double Hvb;
        double Svb;
        double SigmaD;
        double SigmaP;
        double SigmaH;
        double Sigma;
        double Vmol;    //Volume molar @298K [kmol/m³] calculado por MG
        double w;       //Fator acêntrico calculado por MG

        //insaturados PO
        std::vector<std::vector<int> > mGCH2dCH;
        std::vector<std::vector<int> > mGCH2dC;
        std::vector<std::vector<int> > mGCHdCH;
        std::vector<std::vector<int> > mGCdC;
        std::vector<std::vector<int> > mGCHdC;
        std::vector<std::vector<int> > mGCH2dCdCH;
        std::vector<std::vector<int> > mGCH2dCdC;  
        std::vector<std::vector<int> > mGCHdCdCH;

        //Auxiliares reacionais
        // unsigned int Het;       //Definição de heteroátomos
        // unsigned int Aro;       //Definição de aromaticos
        // unsigned int Naf;       //Defininção de Naftênicos
    };

    const std::vector < double > ParamAjustMG {
                                    144.0977,   // 0 - Tm0 em K
                                    244.7889,   // 1 - Tb0 em K
                                    181.6738,   // 2 - Tc em K 
                                    0.0519,     // 3 - Pc1 em bar 
                                    0.1155,     // 4 - Pc2 em bar^-0,5 
                                    14.6182,    // 5 - Vc0 em cm3/mol 
                                    8.5016,     // 6 - Gf0 em KJ/mol 
                                    83.9657,    // 7 - Hf0 em KJ/mol
                                    10.4327,    // 8 - Hv0 em KJ/mol
                                    -1.2993,    // 9 - Hfus0 em KJ/mol 
                                    0.7520,     // 10 - LogKow0 em [-]
                                    150.0218,   // 11 - Fp0 em k
                                    71.2584,    // 12 - Ait1^b em [-]
                                    525.93,     // 13 - Ait2^b em K
                                    15.0884,    // 14 - Hvb0 em KJ/mol
                                    83.7779,    // 15 - Svb0 em KJ/mol
                                    20.7339,    // 16 - sigma0 em MPa^1/2
                                    0.9132,     // 17 - wA em [-]
                                    0.0447,     // 18 - wB em [-]
                                    1.0039,     // 19 - wC em [-]
                                    0.0123};    // 20 - Vm0 em cm3/Kmol 
    //

    //                                  Tb , Tc, Pc, Vc, Tm, Gf, Hf, Hfus, LogKow, Fp, Hv, Hvb, Svb, SigmaD, SigmaP, SigmaH, Sigma, AiT1, AiT2, w, Vm  (Simultaneos Regression Method)
    const Eigen::MatrixXd ParMGPO  {{0.9218,1.0898,0.0100,63.9854,0.7555,-20.4345,-84.0390,1.1420,0.1152,33.0909,1.6141,2.3797,-0.4840,7.5983,2.3037,2.2105,-1.8029,-0.3516,-57.8605,0.0017,0.0238},        // 00 - CH3
                                    {0.5780,3.4604,0.0101,56.8278,0.2966,8.1877,-20.6506,2.6516,0.4594,11.4107,4.8014,2.3004,0.6423,-0.0023,-0.1664,-0.2150,-0.1323,0.1009,2.6047,0.0019,0.0166},           // 01 - CH2
                                    {-0.1189,4.6659,0.0107,37.2813,-0.5960,33.2787,37.6525,-0.9254,0.4300,-17.7416,5.7553,1.6577,0.6894,-7.5390,-3.3851,-2.6826,1.0139,0.5829,79.2115,0.0029,0.0084},       // 02 - CH
                                    {0.7332,3.7648,0.0046,45.2515,0.6317,19.9241,-0.6732,1.4441,0.1960,18.7143,4.0975,2.5866,0.1971,3.1133,0.9190,0.7694,-0.0782,-0.1131,-5.1387,0.0013,0.0120},            // 03 - aCH
                                    {0.7147,3.6308,0.0068,48.7395,0.5313,8.7292,-29.2709,0.9905,0.1646,17.9616,3.6900,2.3438,-0.0399,2.6693,0.5104,0.6159,-0.3085,0.0198,-3.1112,0.0009,0.0160},            // 04 - nCH2
                                    {0.4137,4.7661,0.0065,38.8249,-0.1298,0.9855,-72.8957,0.8855,0.4600,7.5875,4.8950,2.2834,0.6520,-3.6903,-3.2065,-0.5171,-0.6304,0.0719,15.8622,-0.0078,0.0059},         // 05 - nCH
                                    {1.2531,17.9650,0.0033,4.5530,2.0942,28.3757,20.2301,1.6799,0.2148,18.8824,5.7969,4.7372,2.2507,-3.0510,-1.1299,-2.7569,-0.1111,-34.8062,-54.7826,-0.0014,0.0046},      // 06 - aCfaC
                                    {1.1611,18.2562,-0.0068,8.4391,1.6734,29.9678,13.1628,1.5330,0.1295,0.5527,5.5248,3.4588,1.5553,-4.4881,-0.3554,-2.6898,-1.2301,23.2758,-2.4889,-0.0004,0.0035},        // 07 - aCfnC
                                    {1.2616,7.9731,0.0153,99.8977,1.1155,15.9981,-33.2700,3.0939,0.5432,43.3846,6.9939,4.1963,-0.8843,3.0811,-1.1846,0.1713,-0.5759,-0.3481,-214.5329,0.0031,0.0302},       // 08 - aC_CH3
                                    {0.8530,11.8624,0.0162,77.4668,-0.1922,43.9629,31.7820,2.3781,0.6971,24.7144,8.9948,4.5211,1.0856,-5.0406,-4.1990,-3.1805,-0.0100,0.1472,-51.0665,0.0042,0.0229},       // 09 - aC_CH2
                                    {0.0274,12.1443,0.0183,68.6439,-0.8845,60.0040,82.2207,2.2627,1.0146,-4.5611,9.0081,0,0,-13.6223,-5.2318,-7.0890,-0.3146,0.5689,19.7360,0.0020,0.0162},                 // 10 - aC_CH
                                    {-0.6053,7.1924,0.0215,59.4496,-1.2958,103.5315,142.8566,-0.9647,1.3668,-9.6315,13.7692,5.6746,3.4909,-21.2544,-8.0888,-7.8370,1.2302,0.6880,-28.5784,0.0015,0.0095},   // 11 - aC_C
                                    {1.0507,9.5985,-0.0057,23.7909,2.0185,84.8886,54.1833,4.5449,-0.3810,49.5376,11.8029,8.0608,9.8797,3.4095,5.9010,3.9266,-0.8197,-0.0832,-62.5438,0.0008,0.0052},        // 12 - aN
                                    {2.5879,16.1803,0.0115,106.0166,2.1851,6.9435,-46.5497,8.3721,0.3193,66.4735,15.1089,9.9537,2.4892,8.5494,3.8374,4.2415,-0.4082,-0.0764,-87.3884,0.0042,0.0341},        // 13 - CH2SH
                                    {2.2640,12.9218,0.0068,116.9153,3.3490,47.0463,-49.1973,12.4026,-1.4619,58.6201,14.0998,9.6773,5.7817,8.1995,5.2101,6.7984,0.3355,0.0522,36.9000,0.0104,0.0262},        // 14 - CH2NH2
                                    {0.8317,21.0902,-0.0001,23.4172,0.9510,36.6009,49.4849,0.2638,0.2487,29.4370,6.3586,0,0,-6.6777,-3.0946,-1.4992,-1.2037,-1.9551,-5.7471,0.0016,0.0071},                 // 15 - aCqq
                                    {1.7834,28.7804,0.0093,0,0.9894,0,0,0,0.6474,0,0,0,0,-3.5646,-2.0986,-1.3574,0,0,0,0.0053,0.0317},                                                                      // 16 - aCS
                                    {1.0680,18.8249,0.0086,96.0454,1.4958,179.8623,170.1188,3.2332,0.3297,25.8875,18.6549,0,0,0,0,0,1.4103,1.0105,132.0401,0.0049,0.0190},                                  // 17 - aCN
                                    {1.1476,13.0925,0.0074,31.4549,1.3231,-70.0222,-65.8355,1.2779,0.0450,2.3310,14.0271,8.5679,6.9671,-4.4171,-1.2975,0.0149,1.8991,1.3612,-25.1637,0.0068,0.0084},        // 18 - aCO
                                    {2.0774,13.1905,-0.0079,68.4257,2.6727,-3.1011,-19.2417,5.8759,0.1163,0,11.0794,8.2684,2.2265,8.7339,2.5118,6.0326,0.8350,0,0,0.0022,0.0161},                           // 19 - SH
                                    {1.9867,9.4904,-0.0109,56.7075,3.8058,23.7639,5.1575,6.3670,-1.0653,79.3176,18.4004,11.6741,9.8812,9.0936,6.3555,9.0778,5.1519,0.2561,-2.1921,0.0066,0.0096},           // 20 - NH2
                                    {2.2476,10.1672,-0.0071,24.4092,3.2424,-161.8492,-213.8185,3.9494,-1.3365,87.6576,23.9705,17.3709,21.1170,8.0503,5.2379,11.8005,3.0524,-0.1516,32.3056,0.0180,0.0042},  // 21 - OH
                                    {1.9665,15.1201,0.0178,122.2872,2.8589,-114.0724,-154.6631,8.3352,-0.2250,63.7628,16.7382,9.9376,4.9255,0.5371,1.2706,-0.0788,0.7468,0.1214,94.5070,0.0080,0.0283},     // 22 - CH2CO
                                    {2.6907,14.1929,0.0189,132.0797,3.2535,-149.7355,-213.0782,6.2525,-0.3980,87.7425,12.8497,8.9922,2.8550,8.1107,6.3823,3.4394,0.4030,0.1335,119.8078,0.0077,0.0347},     // 23 - CH3CO
                                    {1.1925,14.7945,0.0163,94.3910,1.5650,-90.5148,-91.1938,5.9492,0.2663,0,16.1213,9.1129,4.9603,0,0,0,2.5361,0,0,0.0091,0.0199},                                          // 24 - CHCO
                                    {2.2519,29.9384,0.0136,81.0181,2.4974,-79.8262,-85.9307,4.2919,-0.1345,73.3748,22.8491,0,0,-3.8648,0.0858,-0.8306,1.3457,8.4848,260.6371,0.0089,0.0208},                // 25 - aC_CO
                                    {2.1021,11.2208,0.0045,58.6962,2.9059,-120.2465,-165.7752,9.7013,-0.9687,72.5327,12.7980,8.7662,4.9348,7.8411,7.8726,5.3761,1.4631,0.4254,-30.4259,0.0079,0.0167},      // 26 - aldeido_CO
                                    {2.6711,25.4355,0.0138,98.8795,3.1633,-89.2637,-119.3350,8.4880,-0.1943,91.5521,0,0,0,3.7273,3.7488,2.3964,1.2651,6.0646,-33.3858,0.0102,0.0187},                       // 27 - aldeido_aC_CO
                                    {-0.6495,6.6169,0.0075,25.0561,-0.3679,66.0035,96.1819,-2.0279,0.8143,-36.6949,4.9183,0.6634,-0.1297,-15.6455,-5.1979,-6.4821,1.2449,0.6668,95.9781,0.0008,-0.0015},    // 28 - C4rio
                                    {1.4953,5.2031,0.0181,106.7610,1.0430,76.0611,23.7826,1.8927,0.2304,42.4673,-0.4710,3.6704,0.3739,7.7504,3.6752,2.7673,-2.1839,-0.2295,-72.0187,0.0027,0.0333},         // 29 - CH2=CH
                                    {1.2001,8.2552,0.0194,97.7190,0.6600,95.8888,82.4195,5.4232,0.4297,10.1964,-2.4067,3.6489,0.1313,0.4284,3.0492,0.8631,-0.3708,0.0768,-116.5855,0.0046,0.0244},          // 30 - CH=CH
                                    {1.0308,7.3554,0.0176,93.6274,0.3327,91.0854,78.1789,3.1782,0.3122,6.0196,-2.9764,3.2053,0.5815,0.1956,2.3059,1.0623,-0.8300,0.0794,-44.1695,0.0020,0.0230},            // 31 - CH2=C
                                    {0.7646,10.0135,0.0170,87.2150,-0.3944,114.7548,138.6480,4.3549,0.6026,-15.4082,-5.6340,4.7812,1.9897,-7.0086,1.2790,-0.1204,1.6192,0.5499,-66.0852,0.0036,0.0108},     // 32 - CH=C
                                    {0.4080,13.5316,0.0209,101.2116,-0.9826,142.3553,201.5682,2.7754,0.5563,-44.0014,-7.2229,5.0939,3.0012,-14.6160,-1.4590,-2.9995,2.7550,0.8132,11.7928,0.0007,-0.0021},  // 33 - C=C
                                    {2.1764,11.0144,0.0261,140.9994,1.6865,209.8955,161.8935,6.0995,0,0,11.1663,6.5451,1.2959,7.7662,2.2910,3.6733,-2.7895,0,0,0.0037,0.0431},                              // 34 - CH2=C=CH
                                    {1.7104,12.6575,0.0261,148.4107,1.5154,229.8483,213.1922,7.0169,0,0,0,0,0,-0.0965,-2.1075,0.0791,-1.3229,0,0,0.0020,0.0365},                                            // 35 - CH2=C=C
                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},                                                                                                                            // 36 - CH=C=CH
                                    {1.3771,6.4667,0.0091,83.8235,1.2254,75.1097,37.0718,1.0124,0.2761,26.9856,7.4405,4.8318,1.2589,7.1860,0.7482,1.9775,-0.7844,-0.4279,-49.4820,0.0024,0.0246},           // 37 - nCH=CH
                                    {0.9247,9.1574,0.0087,74.9435,0.8270,63.4709,41.6830,3.7399,0.6105,33.0981,5.9475,6.2169,4.4891,-1.1666,3.3074,-1.5025,1.1486,0.2903,33.1340,0.0059,0.0020},            // 38 - nCH=C
                                    {0.6624,0,0,0,1.0683,0,0,0,0.8714,0,0,0,0,-7.5565,-7.2796,-17.0120,-1.9522,0,0,0.0127,-0.0361}};                                                                        // 39 - nC=C
    //


    //Matriz com a correlacao entre os grupos de Joback e os grupos MSOL extendido
    //				  			          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41
    //				 			         A6 A4 A3 A2 N6 N5 4a 4n 3a 3m 3n 2n 2a 2m 2f 1a 1m 1n Rp Rm Rn Ra br b2 Mn Ma Ho Hn Aa Am An NS RS ANaANnNN RN NO RO KO Ni VO
    const Eigen::MatrixXd Mjoback_SOL  {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  0  -CH3nr 	Alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-2,-1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0}, //  1  -CH2-nr	Alif�tico
                                        {0, 0, 0, 0, 6, 5, 4, 2, 3, 2, 1, 0, 2, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0, 0,-1,-2,-1, 0, 0,-2,-1, 0,-1, 0, 0, 0, 0}, //  2  >CH2r 	Naft�nico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  3  >CH-nr 	Alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 2, 0, 1, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  4  >CH-r 	Naft�nico condensado
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  5  >CH-r 	Naft�nico substitu�do
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  6  =CH-nr 	Olefina linear
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, //  7  =CH-r 	Olefina c�clica
                                        {6, 2, 1, 0, 0, 0,-2, 0,-2,-1, 0, 0,-2,-1,-1,-2,-1, 0, 0, 0, 0,-1, 0, 0, 0,-1, 0, 0,-2,-1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0}, //  8  CH 		Arom�tico
                                        {0, 2, 2, 0, 0, 0, 2, 0, 2, 1, 0, 0, 0, 0,-1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //  9  C 		Arom�tico condensado perif�rico
                                        {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // 10  C 		Arom�tico condensado interno
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // 11  =C<r 	Arom�tico substitu�do
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // 12  -S-r 	Enxofre benzoico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // 13  -S-nr 	Enxofre alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}, // 14  =N-r 	Nitr�genio piridinico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, // 15  -NH-r 	Nitrog�nio pirrolico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, // 16  -NH-nr	Nitrog�nio alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, // 17  -O-r 	Oxig�nio fur�nico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, // 18  -O-nr 	Oxig�nio - eter alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, // 19  >C=Onr 	Carbonil alif�tico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, // 20  -Ni- 	N�quel Porfir�nico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, // 21  -VO- 	Van�dio Porfir�nico
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};// 22  >C<      Alifático quaternário
    //

    //Matriz com a contribuição Joback para as propriedades
    //				  			                   0     1       2     3     4  
    //				 			                  Tc     Pc      Vc    Tb    Tm
    const Eigen::MatrixXd MContribuicoes_Joback {{0.0141,-0.0012,65   ,23.58,-5.10},    //  0  -CH3nr 	Alif�tico
                                                 {0.0189,0      ,56   ,22.88,11.27},    //  1  -CH2-nr	Alif�tico
                                                 {0.01  ,0.0025 ,48   ,27.15,7.75},     //  2  >CH2r 	Naft�nico
                                                 {0.0164,0.002  ,41   ,21.74,12.64},    //  3  >CH-nr 	Alif�tico
                                                 {0.0122,0.0004 ,38   ,21.78,19.88},    //  4  >CH-r 	Naft�nico condensado
                                                 {0.0122,0.0004 ,38   ,21.78,19.88},    //  5  >CH-r 	Naft�nico substitu�do
                                                 {0.0129,-0.0006,46   ,24.96,8.73},     //  6  =CH-nr 	Olefina linear
                                                 {0.0082,0.0011 ,41   ,26.73,8.13},     //  7  =CH-r 	Olefina c�clica
                                                 {0.0082,0.0011 ,41   ,26.73,8.13},     //  8  CH 		Arom�tico
                                                 {0.0042,0.0061 ,27   ,21.32,60.15},    //  9  C 		Arom�tico condensado perif�rico
                                                 {0.0042,0.0061 ,27   ,21.32,60.15},    // 10  C 		Arom�tico condensado interno
                                                 {0.0042,0.0061 ,27   ,21.32,60.15},    // 11  =C<r 	Arom�tico substitu�do
                                                 {0.0019,0.0051 ,38   ,52.1 ,79.93},    // 12  -S-r 	Enxofre benzoico
                                                 {0.0119,0.0049 ,54   ,68.78,34.4},     // 13  -S-nr 	Enxofre alif�tico
                                                 {0.0085,0.0076 ,34   ,57.55,68.4},     // 14  =N-r 	Nitr�genio piridinico
                                                 {0.013 ,0.0114 ,29   ,52.82,101.51},   // 15  -NH-r 	Nitrog�nio pirrolico
                                                 {0.0295,0.0077 ,35   ,50.17,52.66},    // 16  -NH-nr	Nitrog�nio alif�tico
                                                 {0.0098,0.0048 ,13   ,31.22,23.05},    // 17  -O-r 	Oxig�nio fur�nico
                                                 {0.0168,0.0015 ,18   ,22.42,22.23},    // 18  -O-nr 	Oxig�nio - eter alif�tico
                                                 {0.038 ,0.0031 ,62   ,76.75,61.2},     // 19  >C=Onr 	Carbonil alif�tico
                                                 {0     ,0      ,0    ,0    ,0},        // 20  -Ni- 	N�quel Porfir�nico
                                                 {0     ,0      ,0    ,0    ,0},        // 21  -VO- 	Van�dio Porfir�nico
                                                 {0.0067,0.0043 ,27   ,18.25,46.43}};   // 22  >C<      Alifático quaternário
    //

    //Matriz com o numero de átomos por grupo joback
    //				  			           0
    //				 			           NAtomos
    const Eigen::MatrixXd MAtomos_Joback {{4},  //  0  -CH3nr 	Alif�tico
                                          {3},  //  1  -CH2-nr	Alif�tico
                                          {3},  //  2  >CH2r 	Naft�nico
                                          {2},  //  3  >CH-nr 	Alif�tico
                                          {2},  //  4  >CH-r 	Naft�nico condensado
                                          {2},  //  5  >CH-r 	Naft�nico substitu�do
                                          {2},  //  6  =CH-nr 	Olefina linear
                                          {2},  //  7  =CH-r 	Olefina c�clica
                                          {2},  //  8  CH 		Arom�tico
                                          {1},  //  9  C 		Arom�tico condensado perif�rico
                                          {1},  // 10  C 		Arom�tico condensado interno
                                          {1},  // 11  =C<r 	Arom�tico substitu�do
                                          {1},  // 12  -S-r 	Enxofre benzoico
                                          {1},  // 13  -S-nr 	Enxofre alif�tico
                                          {1},  // 14  =N-r 	Nitr�genio piridinico
                                          {2},  // 15  -NH-r 	Nitrog�nio pirrolico
                                          {2},  // 16  -NH-nr	Nitrog�nio alif�tico
                                          {1},  // 17  -O-r 	Oxig�nio fur�nico
                                          {1},  // 18  -O-nr 	Oxig�nio - eter alif�tico
                                          {2},  // 19  >C=Onr 	Carbonil alif�tico
                                          {1},  // 20  -Ni- 	N�quel Porfir�nico
                                          {2},  // 21  -VO- 	Van�dio Porfir�nico
                                          {1}}; // 22  >C<      Alifático quaternário
    //


    // Cosntruindo vetor com nome das propriedades
    const std::vector < std::string > NomesProps
    {
        "MW",                   //  0 - Massa molar mistura
        "densidade",            //  1 - Massa específica da mistura
        "C_hsno",               //  2 - CHSNO: Percentual mássico de Carbono - Análise Elementar
        "c_H_sno",              //  3 - CHSNO: Percentual mássico de Hidrogênio - Análise Elementar
        "ch_S_no",              //  4 - CHSNO: Percentual mássico de Enxofre - Análise Elementar
        "chs_N_o",              //  5 - CHSNO: Percentual mássico de Nitrogênio - Análise Elementar
        "chsn_O",               //  6 - CHSNO: Percentual mássico de Oxigênio - Análise Elementar
        "S_ara",                //  7 - SARA: Percentual mássico de Saturados - SARA
        "s_A_ra",               //  8 - SARA: Percentual mássico de Aromáticos - SARA
        "sa_R_a",               //  9 - SARA: Percentual mássico de Resinas - SARA
        "sar_A",                // 10 - SARA: Percentual mássico de Asfaltenos - SARA
        "RMN_Cinsat",           // 11 - RMNC: Percentual molar de Carbono insaturado - RMNC
        "RMN_Csat",             // 12 - RMNC: Percentual molar de Carbono Saturado - RMNC
        "RMN_CaroAlquil",       // 13 - RMNC: Percentual molar de Carbono Aromático substituído por alquila - RMNC
        "RMN_CaroH",            // 14 - RMNC: Percentual molar de Carbono Aromático não substituído, com Hidrogênio - RMNC
        "RMN_fa",               // 15 - RMNC: Fator de aromaticidade - RMNC
        "RMN_Hmar",             // 16 - RMNH: Percentual molar de Hidrogênio em anéis mono aromáticos - RMNH 
        "RMN_Hdar",             // 17 - RMNH: Percentual molar de Hidrogênio em anéis aromáticos condensados (diaromáticos ou maiores) - RMNH
        "RMN_Haro",             // 18 - RMNH: Percentual molar de Hidrogênio aromático total (Hmar + Hdar) - RMNH
        "RMN_Holef",            // 19 - RMNH: Percentual molar de Hidrogênio olefínico - RMNH    
        "RMN_Halpha",           // 20 - RMNH: Percentual molar de Hidrogênio saturado em posição alpha - RMNH
        "RMN_Hbeta1",           // 21 - RMNH: Percentual molar de Hidrogênio saturado em posição beta + CHs de cicloparafinas + CHs de cadeias parafínicas - RMNH
        "RMN_Hbeta2",           // 22 - RMNH: Percentual molar de Hidrogênio saturado em posição maior que beta + CH2 de cadeias parafínicas + CH2 de cicloparafinas - RMNH 
        "RMN_Hbeta",            // 23 - RMNH: Percentual molar de Hidrogênio saturado beta total (Hbeta_1 + Hbeta_2) - RMNH
        "RMN_Hgama",            // 24 - RMNH: Percentual molar de Hidrogênio saturado de grupos CH3 terminais ou isolados
        "RMN_Hsat",             // 25 - RMNH: Percentual molar de Hidrogênio saturado total (Halpha + Hbeta + Hgamma)
        "RMN_RelacaoCH",        // 26 - RMN: Razão molar entre C e H - RMN
        "RMN_AlcLinear",        // 27 - RMN:
        "RMN_TMC",              // 28 - RMN: Tamanho médio de Cadeias alifáticas - RMN
        "RMN_CH3terminal",      // 29 - RMN: Percentual molar de Carbono metílico (CH3) terminal
        "RMN_CH3ramif",         // 30 - RMN: Percentual molar de Carbono métílico (CH3) em ramificação
        "RMN_CH3Tot_CH2Par",    // 31 - RMN: Relação molar entre CH3Total e CH2parafínico
        "RMN_CH3Ram_CH2Par",    // 32 - RMN: Relação molar entre CH3Ramificado e CH2parafínico
        "BPM",                  // 33 - CAP: Percentual molar de moléculas de Baixo MW (MW<425)
        "MPM",                  // 34 - CAP: Percentual molar de moléculas de Médio MW (425<MW<3700)
        "APM",                  // 35 - CAP: Percentual molar de moléculas de Alto MW (MW>3700)
        "GCMS_nPar",            // 36 - GCMS:
        "GCMS_isoPar",          // 37 - GCMS:
        "GCMS_Par",             // 38 - GCMS:
        "GCMS_cPar",            // 39 - GCMS:  
        "GCMS_MonoAro",         // 40 - GCMS:
        "GCMS_DiAro",		    // 41 - GCMS:
        "GCMS_TriAro",	        // 42 - GCMS:
        "GCMS_TetraAro",        // 43 - GCMS:
        "GCMS_PentaAro",	    // 44 - GCMS:
        "GCMS_Sulfur",          // 45 - GCMS:
        "SFC_Sat",              // 46 - SFC: 
        "SFC_MonoAro",          // 47 - SFC: 
        "SFC_DiAro",            // 48 - SFC:
        "SFC_PoliAro",          // 49 - SFC:
        "SFC_TotalAro",         // 50 - SFC:
        "MCR",                  // 51 - MCR:
        "wtASPH",               // 52 - Teor de asfalteno - insoluveis em nC7
        "d156",                 // 53 - d15,6°C/15,6°C
        "API",                  // 54 - °API
        "Visc378",              // 55 - Viscosidade cinemática @ 37,8°C mm²/s = cSt
        "Visc989",              // 56 - Viscosidade cinemática @ 98,9°C mm²/s = cSt
        "Avisc",                // 57 - Parâmetro A da equação de Walther-ASTM para viscosidade
        "Bvisc",                // 58 - Parâmetro B da equação de Walther-ASTM para viscosidade
        "IR",                   // 59 - Índice de refração @20°C
        "PtoAnilina",           // 60 - Ponto de anilina K
        "PtoFuligem",           // 61 - Ponto de fuligem mm
        "PtoConge",             // 62 - Ponto de Congelamento K
        "PtoNevoa",             // 63 - Ponto de Nevoa K
        "PtoFluidez",           // 64 - Ponto de Fluidez K
        "Visc60",               // 65 - Viscosidade cinemática @ 60,0°C mm²/s = cSt - 55
        "Visc822",              // 66 - Viscosidade cinemática @ 82,2°C mm²/s = cSt - 55
        "Visc100",              // 67 - Viscosidade cinemática @ 100,0°C mm²/s = cSt - 55
        "Visc135",              // 68 - Viscosidade cinemática @ 135,0°C mm²/s = cSt - 55
        "Visc150",              // 69 - Viscosidade cinemática @ 150,0°C mm²/s = cSt - 55
        "CoefTeb",
        "PEVm_T0_K",
        "PEVm_T5_K",
        "PEVm_T10_K",
        "PEVm_T20_K",
        "PEVm_T30_K",
        "PEVm_T40_K",
        "PEVm_T50_K",
        "PEVm_T60_K",
        "PEVm_T70_K",
        "PEVm_T80_K",
        "PEVm_T90_K",
        "PEVm_T95_K",
        "PEVm_T100_K"
    };


    //Definicao das dimensões das variáveis
    const unsigned int MAGcolumns = 5;          //numero de colunas da MAG
    const unsigned int MAGmaxline = 600;        //número de linhas da MAG
    const unsigned int maxLigas_e = 4;          //número máximo de ligações que um elemento pode fazer na molécula
    const unsigned int maxNucleos = 5;          //Número máximo de nucleo por molecula
    const unsigned int nParametros = 14;        //numero de parâmnetros RE
    const unsigned int numercasasdecimais = 4;  //número de casas decimais para salvar arquivos de prorpruiedades e parâmetros
    const unsigned int Nucleosize = maxNucleos + 1;      //Número máximo de nucleo por molecula
    const unsigned int NligasMaxPorAtomo = maxNucleos - 1;
    const unsigned int SOLexpsize = Mjoback_SOL.cols() + 2;
    const unsigned int MGPOsize = ParMGPO.rows();

    const Eigen::MatrixXd T_Mjoback_SOL = Mjoback_SOL.transpose();
};

#endif // CONSTANTES_HPP