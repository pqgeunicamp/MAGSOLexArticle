#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <vector>
#include <cmath> //divisão e arredondamento usada em craqueamento
#include <sstream>
#include <random>
#include <algorithm>
#include <string>

// reconhecimento da classe dos seguintes arquivos
#include "./repomme/Mistura.hpp"
#include "./repomme/MG.hpp"
#include "./repomme/Randomic.hpp"

int main(int argc, char** argv)
{
    unsigned int user_option = 0;
    // Lendo argumentos de entrada
	for (int i = 1 ; i < argc ; i++)
	{
		// Se Exemplo 1 - calculo de propriedades
		if(strcmp(argv[i], "-Properties") == 0) user_option = 1;
        else if(strcmp(argv[i], "-Reactions") == 0) user_option = 2;
	}

    // Informações iniciais
    std::cout << std::endl;
    std::cout << "MAG-SOLex molecular representation: a methodology for handling complex molecules in algorithms" << std::endl;
    std::cout << "Example of using the MAG-SOLex representation." << std::endl;
    std::cout << "by Diego Telles, Karina Costa & Tarcisio Dantas" << std::endl;
    std::cout << "PETROBRAS - CENPES" << std::endl;
    std::cout << "UNICAMP - DEPro" << std::endl;
    std::cout << __DATE__ << " as " << __TIME__ << std::endl;
    std::cout << std::endl;

    // Cria objeto que contem as informações sobre a mistura de moléculas
    Mistura MMix;

    if (user_option == 1)
    {// Properties 

        #pragma region Introduction

        std::cout << "Properties calculation example" << std::endl;
        std::cout << std::endl;
    
        MMix.CarregaMistura("./run/ExProperties", "Composicao.csv", "SOLex.csv", "MAG.csv");


        bool printMolecules = true;
        if (printMolecules)
        {
            std::cout << "Introducing the molecules of the properties calculation example using the MAG and SOLex representations" << std::endl;
            std::cout << std::endl;
            MMix.ImprimeMAG();
            MMix.ImprimeSOLex();
        }
        else
        {
            std::cout << "If you want to see molecules set true in printMolecules variable" << std::endl;
            std::cout << std::endl;
        }

        std::cout << "SOLex matrix in eigenformat - molecules in the rows [SOLexEigen]:" << std::endl;
        std::cout << "A6,A4,A3,A2,N6,N5,N4a,N4n,N3a,N3m,N3n,N2n,N2a,N2m,N2f,N1a,N1m,N1n,Rp,Rm,Rn,Ra,br,br2,MEn,MEa,IHo,IHn,AAa,AAm,AAn,NS,RS,ANa,ANn,NN,RN,NO,RO,KO,Ni,VO" << std::endl;
        std::cout << MMix.SOLexEigen << std::endl;
        std::cout << std::endl;

        std::cout << "Matrix relating Joback and Reid groups to SOLex attributes - Joback groups in the rows [Mjoback_SOL] (Supporting Information - Table 1):" << std::endl;
        std::cout << "A6,A4,A3,A2,N6,N5,N4a,N4n,N3a,N3m,N3n,N2n,N2a,N2m,N2f,N1a,N1m,N1n,Rp,Rm,Rn,Ra,br,br2,MEn,MEa,IHo,IHn,AAa,AAm,AAn,NS,RS,ANa,ANn,NN,RN,NO,RO,KO,Ni,VO" << std::endl;
        std::cout << CTE::Mjoback_SOL << std::endl;
        std::cout << std::endl;

        #pragma endregion Introduction

        #pragma region Joback
        /*
        ######################################################################################
        ####    Calculating Properties Using Joback and Reid Group Contribution Method    ####
        ######################################################################################
        */

        // Calculating Joback groups by molecules 
        MMix.Mjoback_moleculas = MMix.SOLexEigen*CTE::T_Mjoback_SOL;

        std::cout << "Matrix of Joback and Reid groups by molecules - molecules in the rows [Mjoback_moleculas = SOLexEigen*Transpose(Mjoback_SOL)]:" << std::endl; 
        std::cout << "CH3-nr\t-CH2-nr\t>CH2-r\t>CH-nr\tCH-r1\tCH-r2\t=CH-nr\t=CH-r\t=CH-r\t=C<r1\t=C<r2\t=C<r\t-S-r\t=-S-nr\t=N-r\t-NH-r\t-NH-nr\t-O-r\t-O-nr\tC=Onr\t-Ni-\t-VO-\t>C<" << std::endl;
        std::cout << MMix.Mjoback_moleculas << std::endl;
        std::cout << std::endl;

        std::cout << "Matrix with the values of the Joback group contributions for each property - Joback groups in the rows [MContribuicoes_Joback] (JOBACK & REID,1987 - doi.org/10.1080/00986448708960487):" << std::endl; 
        std::cout << "Tc\tPc\tVc\tTb\tTm" << std::endl;
        std::cout << CTE::MContribuicoes_Joback << std::endl;
        std::cout << std::endl;

        // calculating sum of contributions Joback per molecule
        Eigen::MatrixXd MContJoback_mol = MMix.Mjoback_moleculas*CTE::MContribuicoes_Joback;

        std::cout << "Matrix with the sum of contributions per molecule for each property - molecules in the rows [MContJoback_mol = Mjoback_moleculas*MContribuicoes_Joback]:" << std::endl; 
        std::cout << "Tc\tPc\tVc\tTb\tTm" << std::endl;
        std::cout << MContJoback_mol << std::endl;
        std::cout << std::endl;

        // calculating the number of atoms per molecule
        Eigen::MatrixXd MNa = MMix.Mjoback_moleculas*CTE::MAtomos_Joback;
        // std::cout << "MNa:" << std::endl;
        // std::cout << MNa << std::endl;
        // std::cout << std::endl;

        // calculating the property value per molecule
        Eigen::MatrixXd JobackProperties_mol(MMix.nMol,MContJoback_mol.cols());

        for (unsigned int i = 0 ; i < JobackProperties_mol.rows(); i++)
        {
            JobackProperties_mol(i,3) = 198.2 + MContJoback_mol(i,3);                                                                           //Tb [K]
            JobackProperties_mol(i,0) = JobackProperties_mol(i,3)*pow((0.584 + 0.965*MContJoback_mol(i,0) - pow(MContJoback_mol(i,0),2)),-1);   //Tc [K]
            JobackProperties_mol(i,1) = pow((0.113 + 0.0032*MNa(i,0) - MContJoback_mol(i,1)),-2);                                                     //Pc [bar]
            JobackProperties_mol(i,2) = 17.5 + MContJoback_mol(i,2);                                                                            //Vc [cm3/mol]
            JobackProperties_mol(i,4) = 122.5 + MContJoback_mol(i,4);                                                                           //Tm [K]
        }

        std::cout << "Matrix with the calculated Joback properties per molecule - molecules in the rows [JobackProperties_mol]:" << std::endl; 
        std::cout << "Tc\tPc\tVc\tTb\tTm" << std::endl;
        std::cout << JobackProperties_mol << std::endl;
        std::cout << std::endl;
        
        #pragma endregion Joback

        #pragma region Marrero
        /*
        #######################################################################################
        ####    Calculating Properties Using Marrero and Gani Group Contribution Method    ####
        #######################################################################################
        */

        /*
        Calculating the quantity of first-order Marrero and Gani groups for the molecules. 
        In this section, an algorithm similar to that presented in Algorithm 1 of the article 
        is applied, where the MAG and SOLex representations are used together. The information 
        contained in SOLex is used to guide the inspection functions in MAG.
        */

        bool controle = 0;
        for (unsigned int m = 0 ; m < MMix.mol.size() ; m++)
        {

            bool Aro = false;
            bool Naf = false;
            bool Olef = false;
            bool Oxi = false;
            bool Nit = false;
            bool Sul = false;

            if (MMix.mol[m].SOLex[0][0] + MMix.mol[m].SOLex[0][1] + MMix.mol[m].SOLex[0][2] + 
                MMix.mol[m].SOLex[0][3] > 0 ) Aro = true;
            if (MMix.mol[m].SOLex[0][4] + MMix.mol[m].SOLex[0][5] + MMix.mol[m].SOLex[0][6] +
                MMix.mol[m].SOLex[0][7]  + MMix.mol[m].SOLex[0][8]  + MMix.mol[m].SOLex[0][9] +
                MMix.mol[m].SOLex[0][10] + MMix.mol[m].SOLex[0][11] + MMix.mol[m].SOLex[0][12] +
                MMix.mol[m].SOLex[0][13] + MMix.mol[m].SOLex[0][14] + MMix.mol[m].SOLex[0][15] +
                MMix.mol[m].SOLex[0][16] + MMix.mol[m].SOLex[0][17] > 0 ) Naf = true;
            if (MMix.mol[m].SOLex[0][26] + MMix.mol[m].SOLex[0][27] != 0) Olef = true;
            if (MMix.mol[m].SOLex[0][37] + MMix.mol[m].SOLex[0][38] + MMix.mol[m].SOLex[0][39] +
                MMix.mol[m].SOLex[0][41] > 0) Oxi = true;
            if (MMix.mol[m].SOLex[0][33] + MMix.mol[m].SOLex[0][34] + MMix.mol[m].SOLex[0][35] +
                MMix.mol[m].SOLex[0][36] > 0) Nit = true;
            if (MMix.mol[m].SOLex[0][31] + MMix.mol[m].SOLex[0][32] > 0) Sul = true;

            
            MMix.mol[m].mGCH2dCH = MG::GCH2dCH(m, MMix);
            MMix.mol[m].mGCH2dC = MG::GCH2dC(m, MMix);
            MMix.mol[m].mGCHdCH = MG::GCHdCH(m, MMix);
            MMix.mol[m].mGCHdC = MG::GCHdC(m, MMix);
            MMix.mol[m].mGCdC = MG::GCdC(m, MMix);
            MMix.mol[m].mGCH2dCdCH = MG::GCH2dCdCH(m, MMix);
            MMix.mol[m].mGCH2dCdC = MG::GCH2dCdC(m, MMix);
            MMix.mol[m].mGCHdCdCH = MG::GCHdCdCH(m, MMix);

            // Clearing the matrix of the quantity of Marrero groups per molecule
            for (unsigned int i = 0 ; i < MMix.MGgroupsPO.cols() ; i++) { MMix.MGgroupsPO(m,i) = 0; }

            // Identifying the groups atom by atom
            // The MAGnode vector indicates the number of atoms, excluding hydrogen, to which atom A is bonded.
            for (unsigned int A = 0 ; A < MMix.mol[m].MAGnode.size() ; A++)
            {
                controle = 0;
                if (MMix.mol[m].MAGnode[A] == 1)
                {
                    if (MG::CH3(m, A, MMix) == 1) MMix.MGgroupsPO(m,0) += 1;
                    else controle = 1;
                }
                else if (MMix.mol[m].MAGnode[A] == 2)
                {
                    if (MG::CH2(m, A, MMix) == 1) MMix.MGgroupsPO(m,1) += 1;
                    else if (MG::nCH2(m, A, MMix) == 1) MMix.MGgroupsPO(m,4) += 1;
                    else if (Aro && MG::aCH(m, A, MMix) == 1) MMix.MGgroupsPO(m,3) += 1;
                    else if (Olef && MG::CH2dCH(m,A,MMix) == 1) MMix.MGgroupsPO(m,29) +=1;
                    else if (Olef && MG::CHdCH(m, A, MMix) == 1) MMix.MGgroupsPO(m,30) +=1;
                    else if (Olef && MG::CH2dCdCH(m, A, MMix) == 1) MMix.MGgroupsPO(m,34) +=1;
                    else if (Olef && MG::CH2dCdC(m, A, MMix) == 1) MMix.MGgroupsPO(m,35) +=1;
                    else if (Olef && MG::CHdCdCH(m, A, MMix) == 1) MMix.MGgroupsPO(m,36) +=1;
                    
                    else controle = 1;
                }
                else if (MMix.mol[m].MAGnode[A] == 3)
                {
                    if (MG::CH(m, A, MMix) == 1) MMix.MGgroupsPO(m,2) += 1;
                    if (MG::nCH(m, A, MMix) == 1) MMix.MGgroupsPO(m,5) += 1;
                    else if (Aro && MG::aCfaC(m, A, MMix) == 1) MMix.MGgroupsPO(m,6) += 1;
                    else if (Aro && MG::aC_CH3(m, A, MMix) == 1) MMix.MGgroupsPO(m,8) += 1;
                    else if (Aro && MG::aC_CH2(m, A, MMix) == 1) MMix.MGgroupsPO(m,9) += 1;
                    else if (Aro && MG::aC_CH(m, A, MMix) == 1) MMix.MGgroupsPO(m,10) += 1;
                    else if (Aro && MG::aC_C(m, A, MMix) == 1) MMix.MGgroupsPO(m,11) += 1;
                    else if (Aro && MG::aCqq(m, A, MMix) == 1) MMix.MGgroupsPO(m,15) += 1;
                    else if (Aro && Naf && MG::aCfnC(m, A, MMix) == 1) MMix.MGgroupsPO(m,7) += 1;
                    else if (Aro && Sul && MG::aCS(m, A, MMix) == 1) MMix.MGgroupsPO(m,16) += 1;
                    else if (Aro && Nit && MG::aCN(m, A, MMix) == 1) MMix.MGgroupsPO(m,17) += 1;
                    else if (Aro && Oxi && MG::aCO(m, A, MMix) == 1) MMix.MGgroupsPO(m,18) += 1;
                    else if (Olef && MG::CH2dC(m, A, MMix) == 1) MMix.MGgroupsPO(m,31) +=1;
                    else if (Olef && MG::CHdC(m, A, MMix) == 1) MMix.MGgroupsPO(m,32) +=1;
                    else if (Olef && MG::CdC(m, A, MMix) == 1) MMix.MGgroupsPO(m,33) +=1;
                    else controle = 1;
                    
                }
                else if (MMix.mol[m].MAGnode[A] == 4)
                {
                    if (MG::C4rio(m, A, MMix) == 1) MMix.MGgroupsPO(m,28) += 1;
                    else controle = 1;
                }
                else if (Sul && MMix.mol[m].MAGnode[A] == 1000)
                {
                    if (MG::CH2SH(m, A, MMix) == 1) MMix.MGgroupsPO(m,13) += 1;
                    else if (MG::SH(m, A, MMix) == 1) MMix.MGgroupsPO(m,19) += 1;
                    else controle = 1;
                }
                else if (Nit && MMix.mol[m].MAGnode[A] == 1001)
                {
                    if (MG::aN(m, A, MMix) == 1) MMix.MGgroupsPO(m,12) += 1;
                    else if (MG::CH2NH2(m, A, MMix) == 1) MMix.MGgroupsPO(m,14) += 1;
                    else if (MG::NH2(m, A, MMix) == 1) MMix.MGgroupsPO(m,20) += 1;
                    else controle = 1;
                }
                else if (Oxi && MMix.mol[m].MAGnode[A] == 1002)
                {
                    if (MG::OH(m, A, MMix) == 1) MMix.MGgroupsPO(m,21) += 1;
                    else controle = 1;
                }
                else if (Oxi && MMix.mol[m].MAGnode[A] == 1003)
                {
                    if (MG::aldeido_CO(m, A, MMix) == 1) MMix.MGgroupsPO(m,26) += 1;
                    else controle = 1;
                }
                
                if (Oxi && controle == 1)
                {
                    if (MG::CH2CO(m, A, MMix) == 1) MMix.MGgroupsPO(m,22) += 1;
                    else if (MG::CH3CO(m, A, MMix) == 1) MMix.MGgroupsPO(m,23) += 1;
                    else if (MG::CHCO(m, A, MMix) == 1) MMix.MGgroupsPO(m,24) += 1;
                    else if (MG::aC_CO(m, A, MMix)  == 1) MMix.MGgroupsPO(m,25) += 1;
                    else if (MG::aldeido_aC_CO(m, A, MMix) == 1) MMix.MGgroupsPO(m,27) += 1;
                }
            }

            int n = 1;
            int nContabilizados = 0;
            for (unsigned int i = 0 ; i < MMix.MGgroupsPO.cols() ; i++)
            {
                if (((i >= 8) && (i <= 11)) || (i == 13) || (i == 14) || (i == 16) || (i == 17) || (i == 18) || (i == 26) || ((i >= 29) && (i <= 33))) n = 2;
                else if ((i == 22) || (i == 23) || (i == 24) || (i == 25) || (i == 27) || ((i >= 34) && (i <= 36))) n =3;
                else n = 1;
                nContabilizados += n*MMix.MGgroupsPO(m,i);
            }

            // mostrando erro de contabilização dos elementos
            if ((nContabilizados != (MMix.mol[m].UltimoC + 1)) && (MMix.aviso == 1))
            {
                std::cout << "MG::ContaMG - erro de contabilizacao molecula " << m  << std::endl;
                std::cout << "  molecula " << m << std::endl;
                std::cout << "  nContabilizados != (MMix.mol[m].UltimoC + 1):" << nContabilizados << "_vs_" << (MMix.mol[m].UltimoC + 1) << m << std::endl;
                
                for (unsigned int i = 0 ; i < MMix.mol[m].MAG.size() ; i++)
                {//varre ligações
                    std::cout  << i << ",";
                    for (unsigned int k = 0 ; k < MMix.mol[m].MAG[i].size(); k++)
                    {//varre colunas
                        std::cout << MMix.mol[m].MAG[i][k];
                        if (k < (MMix.mol[m].MAG[i].size() - 1)) std::cout << ",";
                    }
                    if (MMix.mol[m].MAG[i][0] < 0) i = MMix.mol[m].MAG.size();
                    std::cout << std::endl;
                }
            }
        }

        // Presenting groups counted per molecule.
        std::cout << "Matrix with Marrero and Gani groups counted per molecule - molecules in the rows [MGgroupsPO]:" << std::endl; 
        std::cout << "CH3\tCH2\tCH\taCH\tnCH2\tnCH\taCfaC\taCfnC\taC_CH3\taC_CH2\taC_CH\taC_C\taN\tCH2SH\tCH2NH2\taCqq\taCS\taCN\taCO\tSH\tNH2\tOH\tCH2CO\tCH3CO\tCHCO\taC_CO\taldeido_CO\taldeido_aC_CO\tC4rio\tCH2=CH\tCH=CH\tCH2=C\tCH=C\tC=C\tCH2=C=CH\tCH2=C=C\tCH=C=CH" << std::endl;
        std::cout << MMix.MGgroupsPO << std::endl;
        std::cout << std::endl;

        std::cout << "Matrix with the values of the Marrero and Gani group contributions for each property - Marrero and Gani groups in the rows [ParMGPO] (Hukkerikar et al.,2012 - Supplementary material, Table S5, doi:10.1016/j.ﬂuid.2012.02.010):" << std::endl; 
        std::cout << "Tb\tTc\tPc\tVc\tTm\tGf\tHf\tHfus\tLogKow\tFp\tHv\tHvb\tSvb\tSigmaD\tSigmaP\tSigmaH\tSigma\tAiT1\tAiT2\tw\tVm" << std::endl;
        std::cout << CTE::ParMGPO << std::endl;
        std::cout << std::endl;

        // calculating sum of contributions Marrero per molecule
        Eigen::MatrixXd MContMarrero_mol = MMix.MGgroupsPO*CTE::ParMGPO;

        std::cout << "Matrix with the sum of contributions per molecule for each property - molecules in the rows [MContMarrero_mol = MGgroupsPO*ParMGPO]:" << std::endl; 
        std::cout << "Tb\tTc\tPc\tVc\tTm\tGf\tHf\tHfus\tLogKow\tFp\tHv\tHvb\tSvb\tSigmaD\tSigmaP\tSigmaH\tSigma\tAiT1\tAiT2\tw\tVm" << std::endl;
        std::cout << MContMarrero_mol << std::endl;
        std::cout << std::endl;

        // calculating the property value per molecule
        Eigen::MatrixXd MGProperties_mol(MMix.nMol,MContMarrero_mol.cols());

        for(unsigned int i = 0; i < MMix.mol.size(); i++)
        {
            MGProperties_mol(i,0) = log(MContMarrero_mol(i,0))*CTE::ParamAjustMG[1]; // Normal boiling point (Tb [K])
            MGProperties_mol(i,1) = log(MContMarrero_mol(i,1))*CTE::ParamAjustMG[2]; //Critical temperature (Tc [K])
            MGProperties_mol(i,2) = pow(MContMarrero_mol(i,2) + CTE::ParamAjustMG[4], -2) + CTE::ParamAjustMG[3]; //Critical pressure (Pc [bar])
            MGProperties_mol(i,3) = MContMarrero_mol(i,3) + CTE::ParamAjustMG[5];  //Critical Volum (Vc [cm³/mol])
            MGProperties_mol(i,4) = log(MContMarrero_mol(i,4))*CTE::ParamAjustMG[0]; //Normal melting point (Tm [K])
            MGProperties_mol(i,5) = MContMarrero_mol(i,5) + CTE::ParamAjustMG[6];  //Standard Gibbs energy at 298 K (Gf [kJ/mol])
            MGProperties_mol(i,6) = MContMarrero_mol(i,6) + CTE::ParamAjustMG[7];  //Standard enthalpy of formation at 298 K (Hf [kJ/mol])
            MGProperties_mol(i,7) = MContMarrero_mol(i,7) + CTE::ParamAjustMG[9];  //Standard enthalpy of fusion (Hfus [kJ/mol])
            MGProperties_mol(i,8) = MContMarrero_mol(i,8) + CTE::ParamAjustMG[10]; //Octanol/water particion coefficient
            MGProperties_mol(i,9) = MContMarrero_mol(i,9) + CTE::ParamAjustMG[11];   // Flash point [K]
            MGProperties_mol(i,10) = MContMarrero_mol(i,10) + CTE::ParamAjustMG[8];  //Standard enthalpy of vaporization at 298 K (Hv [kJ/mol])
            MGProperties_mol(i,11) = MContMarrero_mol(i,11) + CTE::ParamAjustMG[14]; // Enthalpy of vaporization at Tb (Hvb [kJ/mol])
            MGProperties_mol(i,12) = MContMarrero_mol(i,12) + CTE::ParamAjustMG[15]; // Entropy of vaporization at Tb (Svb [J/(mol.K)])
            MGProperties_mol(i,13) = MContMarrero_mol(i,13) ; // Hansen Solubility parameter dispersion [MPa^(1/2)]
            MGProperties_mol(i,14) = MContMarrero_mol(i,14) ; // Hansen Solubility parameter polar [MPa^(1/2)]
            MGProperties_mol(i,15) = MContMarrero_mol(i,15) ; // Hansen Solubility parameter h2-bond [MPa^(1/2)]
            MGProperties_mol(i,16) = MContMarrero_mol(i,16) + CTE::ParamAjustMG[16]; // Hildebrand Solubility parameter [MPa^(1/2)]
            MGProperties_mol(i,17) = CTE::ParamAjustMG[12]*pow(10, MContMarrero_mol(i,17)) + CTE::ParamAjustMG[13] + MContMarrero_mol(i,17); // Auto ignition temperature 1 [K]
            MGProperties_mol(i,18) = CTE::ParamAjustMG[12]*pow(10, MContMarrero_mol(i,18)) + CTE::ParamAjustMG[13] + MContMarrero_mol(i,18); // Auto ignition temperature 2 [K]
            MGProperties_mol(i,19) = CTE::ParamAjustMG[17]*(log(pow((MContMarrero_mol(i,19) + CTE::ParamAjustMG[19]),(1/CTE::ParamAjustMG[18])))); // Acentric factor
            MGProperties_mol(i,20) = MContMarrero_mol(i,20) + CTE::ParamAjustMG[20]; //liquid mlar volume at 298 K [cm³/kmol]
        }

        std::cout << "Matrix with the calculated Marrrero and Gani properties per molecule - molecules in the rows [MGProperties_mol]:" << std::endl; 
        std::cout << "Tb\tTc\tPc\tVc\tTm\tGf\tHf\tHfus\tLogKow\tFp\tHv\tHvb\tSvb\tSigmaD\tSigmaP\tSigmaH\tSigma\tAiT1\tAiT2\tw\tVm" << std::endl;
        std::cout << MGProperties_mol << std::endl;
        std::cout << std::endl;

        std::cout << "Figure 9 of the article was constructed using the properties Tb, Tc, Pc and Hf calculated by Marrero e Gani methods for all molecules in this example [MGProperties_mol]." << std::endl; 
        std::cout << std::endl;

        std::cout << "Table 4 of the article was constructed using the quantities of Marrero and Gani groups identified for molecule 5 [MGgroupsPO(5,:)] - 4,6-dimethyl-1-ethylnaphthalene." << std::endl; 
        std::cout << "\tAmount " << std::endl;
        std::cout << "CH3:\t" << MMix.MGgroupsPO(5,0) << std::endl;
        std::cout << "aCH:\t" << MMix.MGgroupsPO(5,3) << std::endl;
        std::cout << "aCfaC:\t" << MMix.MGgroupsPO(5,6) << std::endl;
        std::cout << "aC_CH3:\t" << MMix.MGgroupsPO(5,8) << std::endl;
        std::cout << "aC_CH2:\t" << MMix.MGgroupsPO(5,9) << std::endl;
        std::cout << std::endl; 


        std::cout << "Table 5 of the article was constructed using the properties calculated by both methods for molecule 5 - 4,6-dimethyl-1-ethylnaphthalene." << std::endl; 
        std::cout << "\t\tJoback\tMarrero&Gani " << std::endl;
        std::cout << "Tc[K]:\t\t" << JobackProperties_mol(5,0) << "\t" << MGProperties_mol(5,1) << std::endl;
        std::cout << "Pc[bar]:\t" << JobackProperties_mol(5,1) << "\t" << MGProperties_mol(5,2) << std::endl;
        std::cout << "Vc[cm3/mol]:\t" << JobackProperties_mol(5,2) << "\t" << MGProperties_mol(5,3) << std::endl;
        std::cout << "Tb[K]:\t\t" << JobackProperties_mol(5,3) << "\t" << MGProperties_mol(5,0) << std::endl;
        std::cout << "Tm[K]:\t\t" << JobackProperties_mol(5,4) << "\t" << MGProperties_mol(5,4) << std::endl;
        std::cout << std::endl; 

        #pragma endregion Marrero

        std::cout << "End of the property calculation example" << std::endl;
        std::cout << std::endl; 

    }
    else if (user_option == 2)
    {// Reactions Example
        std::cout << "Reaction calculations example" << std::endl;
        std::cout << std::endl;

        MMix.CarregaMistura("./run/ExReactions", "Composicao.csv", "SOLex.csv", "MAG.csv");

        std::cout << "MAG and SOLex representations of the reactants." << std::endl;
        std::cout << std::endl; 
        MMix.ImprimeMAG();
        MMix.ImprimeSOLex();

        std::cout << "Evaluating the possible reactions:" << std::endl;
        std::cout << std::endl;   

        for (int i = 0; i < MMix.nMol; i++)
        {
            std::cout << "--------------- MOLECULE " << i << " ---------------" << std::endl;
            MMix.CarregaMAGnodes(i);

            #pragma region SatAro
            // --------------------Saturacao de Aromaticos--------------------
            if ((MMix.mol[i].SOLex[0][0] + MMix.mol[i].SOLex[0][1]) > 0) // A6+A4
            {
                // igualando o produto a molecula de origem
                Mistura prod(1);
                prod.mol[0] = MMix.mol[i];

                std::vector<std::vector<int>> AroRings; // identifica os aneis aromaticos
                int A6Ring = -1, A4Ring = -1;
                unsigned int IniciaReac = -1; // Anel A6 ou A4 e primeira ligacao
                int quantasLigacoes = 1;      // Contagem de ligacoes

                for (unsigned int l = 0; l < prod.mol[0].MAG.size(); l++)
                {
                    int NumAnel = prod.mol[0].MAG[l][2];  // numero do primeiro anel
                    int TipoAnel = prod.mol[0].MAG[l][4]; // tipo do primeiro anel

                    // confere se a proxima lig faz parte do mesmo anel, se sim contabiliza, se não push_back em AroRings
                    if (NumAnel == prod.mol[0].MAG[l + 1][2])
                        quantasLigacoes++;
                    else
                    {
                        if (prod.mol[0].MAG[l + 1][2] < 0)
                            break;

                        int InicioAnel = l - quantasLigacoes + 1;                                  // ligacao em que inicia o anel
                        std::vector<int> linha = {NumAnel, quantasLigacoes, TipoAnel, InicioAnel}; // qual anel e quantas linhas
                        AroRings.push_back(linha);
                        quantasLigacoes = 1;
                    }
                }

                // Inicio da reacao
                if (prod.mol[0].SOLex[0][1] > 0) // A4 preferencialmente
                {
                    std::cout << "Applied Reaction: Saturation of A4 Ring: m" << i << " -> m" << MMix.mol.size() << std::endl;
                    // definir qual o anel participa da reação
                    for (unsigned int j = 0; j < AroRings.size(); j++)
                    {
                        if (AroRings[j][1] == 5 && AroRings[j][2] == 2)
                        {
                            A4Ring = AroRings[j][0]; // Seleciona o último A4
                            IniciaReac = AroRings[j][3];
                        }
                    }

                    //*******************************MAG*******************************
                    // alteracao na MAG do produto
                    for (unsigned int k = IniciaReac; k < (IniciaReac + 6); k++)
                        prod.mol[0].MAG[k][4] = 1;

                    //*******************************SOLex*******************************
                    // buscar o anel vizinho de A4Ring para decidir se é N4a ou N4n
                    int Aj_A4Ring = prod.mol[0].MAG[IniciaReac][0];
                    int Viz_A4Ring = -1, Viz_type = -1;

                    for (unsigned int k = 0; k < IniciaReac; k++)
                    {
                        if (prod.mol[0].MAG[k][0] == Aj_A4Ring || prod.mol[0].MAG[k][1] == Aj_A4Ring)
                            Viz_A4Ring = prod.mol[0].MAG[k][2];
                        break;
                    }

                    for (long unsigned int y = 0; y < AroRings.size(); y++)
                    {
                        if (AroRings[y][0] == Viz_A4Ring)
                            Viz_type = AroRings[y][2];
                        break;
                    }

                    prod.mol[0].SOLex[0][1]--; // A4
                    if (Viz_type == 2)
                    {
                        prod.mol[0].SOLex[0][6]++;
                    } // N4a
                    else
                    {
                        prod.mol[0].SOLex[0][7]++;
                    } // N4n

                    for (unsigned int col = 0; col < prod.mol[0].SOLex[1].size() - 2; col++)
                    {
                        prod.mol[0].SOLex[1][col] = prod.mol[0].SOLex[0][col];
                    }
                }

                else if (prod.mol[0].SOLex[0][0] > 0 && prod.mol[0].SOLex[0][1] == 0) // A6 sem A4
                {
                    std::cout << "Applied Reaction: Saturation of A6 Ring: m" << i << " -> m" << MMix.mol.size() << std::endl;

                    //*******************************MAG*******************************
                    for (unsigned int j = 0; j < AroRings.size(); j++)
                    {
                        if (AroRings[j][1] == 6 && AroRings[j][2] == 2)
                        {
                            A6Ring = AroRings[j][0];
                            IniciaReac = AroRings[j][3];
                        }

                        // alteracao na MAG do produto
                        for (unsigned int k = IniciaReac; k < (IniciaReac + 7); k++)
                            prod.mol[0].MAG[k][4] = 1;
                    }

                    //*******************************SOLex*******************************
                    prod.mol[0].SOLex[0][0]--; // A6
                    prod.mol[0].SOLex[0][4]++; // N6
                    if (prod.mol[0].SOLex[0][6] > 0)
                    {
                        prod.mol[0].SOLex[0][6]--; // N4a
                        prod.mol[0].SOLex[0][7]++; // N4n
                    }

                    for (unsigned int col = 0; col < prod.mol[0].SOLex[1].size() - 2; col++)
                    {
                        prod.mol[0].SOLex[1][col] = prod.mol[0].SOLex[0][col];
                    }
                }
                MMix.mol.push_back(prod.mol[0]);

                //If the user wishes the reaction product to continue reacting, they should keep i++ at this point. 
                //Otherwise, remove this command, and the original reagent will proceed to the evaluation of the next reaction rule.
                i++;

                MMix.nMol++;
            }
            #pragma endregion SatAro

            #pragma region Desalquilacao
            // --------------------Dealkylation--------------------
            int Naf = 0, Aro = 0, Cadeia = 0;
            for (unsigned int r = 0; r < 4; r++)
                Aro += MMix.mol[i].SOLex[0][r];
            for (unsigned int r = 4; r < 18; r++)
                Naf += MMix.mol[i].SOLex[0][r];
            Cadeia = MMix.mol[i].SOLex[0][18] + MMix.mol[i].SOLex[0][19] - MMix.mol[i].SOLex[0][24] - MMix.mol[i].SOLex[0][25];

            if (Aro + Naf > 0 && Cadeia > 3)
            {
                std::cout << "Applied Reaction: Dealkylation: m" << i << " -> m" << MMix.mol.size() << " + m" << MMix.mol.size()+1 << std::endl;
                // igualando os produtos a molecula de origem
                Mistura prod(2);
                prod.mol[0] = MMix.mol[i];

                // preenchendo a matriz MAG do produto B com valores de -1
                for (unsigned int l = 0; l < prod.mol[1].MAG.size(); l++)
                {
                    for (unsigned int c = 0; c < 5; c++)
                    {
                        prod.mol[1].MAG[l][c] = -1;
                    }
                }

                //*******************************MAG*******************************
                int C_Ram = 0, iniciaRam = -1, chainSize = 0;
                for (long unsigned int j = 0; j < prod.mol[0].MAG.size(); j++)
                {
                    if (prod.mol[0].MAG[j][0] < 0)
                        break;

                    else if (prod.mol[0].MAG[j][2] == 999)
                        C_Ram++;
                    if (prod.mol[0].MAG[j][2] == 999 && prod.mol[0].MAG[j - 1][2] < 900)
                    {
                        iniciaRam = j;
                    }
                }

                for (int j = iniciaRam + C_Ram; j-- > iniciaRam + 1;)
                {
                    prod.mol[0].MAG.erase(prod.mol[0].MAG.begin() + j); // apaga as linhas da ramificacao
                }

                prod.mol[1].MAG[0][0] = 0; // Ai
                prod.mol[1].MAG[0][1] = 1; // Aj

                unsigned int unsigned_C_Ram = C_Ram;
                for (unsigned int k = 0; k < unsigned_C_Ram - 2; k++)
                {
                    prod.mol[1].MAG[k + 1][0] = prod.mol[1].MAG[k][0] + 1; // Ai
                    prod.mol[1].MAG[k + 1][1] = prod.mol[1].MAG[k][1] + 1; // Aj
                    prod.mol[1].MAG[k][2] = 999;                           // Info
                    prod.mol[1].MAG[k][3] = 0;                             // Core
                    prod.mol[1].MAG[k][4] = 0;                             // Type
                }

                //*******************************SOLex*******************************
                prod.mol[0].SOLex[0][18] = 1;                                    // Rp
                prod.mol[0].SOLex[0][19] = prod.mol[0].SOLex[0][19] - C_Ram + 1; // Rm

                prod.mol[1].SOLex[0][18] = 2;         // Rp
                prod.mol[1].SOLex[0][19] = C_Ram - 1; // Rm

                // corrigindo a linha do core 1
                for (unsigned int col = 0; col < prod.mol[0].SOLex[1].size() - 2; col++)
                {
                    prod.mol[0].SOLex[1][col] = prod.mol[0].SOLex[0][col];
                    prod.mol[1].SOLex[1][col] = prod.mol[1].SOLex[0][col];
                }

                MMix.mol.push_back(prod.mol[0]);
                MMix.nMol++;
                MMix.mol.push_back(prod.mol[1]);
                MMix.nMol++;
            }
            #pragma endregion Desalquilacao
        }

        std::cout << std::endl;
        std::cout << std::endl; 
        std::cout << "MAG and SOLex representations of the set of generated molecules (reactants and products):" << std::endl;
        std::cout << std::endl; 

        std::cout << "Figure 10 of the article" << std::endl;
        MMix.ImprimeMAG();

        std::cout << "Figure 11 of the article" << std::endl;
        MMix.ImprimeSOLex();

        std::cout << "End of the reaction example" << std::endl;
        std::cout << std::endl; 
        
    }
    else 
    {
        std::cout <<"\nExecute the binary by choosing the example to be run.\nFor the propertie calculation example, enter the flag -Properties, \nand for the reaction calculations example, enter the flag -Reactions as follows:\n\n    MAGSOLexExamples -Properties\n\n        or\n\n    MAGSOLexExamples -Reactions\n\nTry again!\n" << std::endl;
    }

    return 0;
}