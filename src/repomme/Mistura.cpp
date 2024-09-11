#include "Mistura.hpp"

//Constructors
Mistura::Mistura()
{
    nMol = -1;
    MW_media.resize(1,1);
}

Mistura::Mistura(int _nMol)
{
    nMol = _nMol;
    mol.resize(nMol);
    z.resize(nMol,1);
    zwt.resize(nMol,1);
    MGgroupsPO.resize(nMol,CTE::MGPOsize);
    MW_media.resize(1,1);

    unsigned int unsigned_nMol = nMol;
    for (unsigned int m = 0 ; m < unsigned_nMol ; m++) 
    {
        mol[m].UltimoC = 0;
        mol[m].SOLex.resize(CTE::Nucleosize);
        for (unsigned int i = 0 ; i < CTE::Nucleosize ; i++) mol[m].SOLex[i].resize(CTE::SOLexpsize);
        mol[m].MAG.resize(CTE::MAGmaxline);
        for (unsigned int i = 0 ; i < CTE::MAGmaxline ; i++) mol[m].MAG[i].resize(CTE::MAGcolumns,-1);
    }
}


void Mistura::CarregaMistura(std::string path, std::string Compfilename, std::string SOLexfilename, std::string MAGfilename)
{

    LeComposicao(path, Compfilename);
    nMol = z.rows();
    mol.resize(nMol);
    LeMAG(path,MAGfilename);
    LeSOLex(path,SOLexfilename);

    CarregaLigacoesPorAtomo();
    CarregaSOLexEigen();

    zwt.resize(nMol,1);
    MGgroupsPO.resize(nMol,CTE::MGPOsize);
    MW_media.resize(1,1);


    bool print = false;
    if (print)
    {
        std::cout << "path:" << path <<std::endl;
        std::cout << "SOLexfilename:" << SOLexfilename <<std::endl;
        std::cout << "MAGfilename:" << MAGfilename <<std::endl;
        std::cout << "Compfilename:" << Compfilename << std::endl << std::endl;

        ImprimeSOLex();
        ImprimeMAG();
        ImprimeComposicao();
        ImprimeMAGnode();
        ImprimeLigacoesPorAtomo();
    }
}

//Lê MAG do arquivo
void Mistura::LeMAG(std::string path, std::string filename)
{
    bool pulaCabecalho = true;

    //std::cout << "Le MAG" << std::endl << std::endl;
    if (nMol < 1)std::cout << "Mistura::LeMAG - Infome o numero de moleculas" << std::endl;
    else
    {
        path += "/" + filename;
        std::ifstream file(path);

        //Msg de Erro de abertura
        if (!file.is_open())
        {
            std::cout << "Mistura::LeMAG - Nao foi possivel abrir " << path << std::endl;
            exit(EXIT_FAILURE);
        }

        //variaveis auxiliares
        std::string line, colum;
        std::vector <int> Ligacao;
        unsigned int unsigned_nMol = nMol;
        Ligacao.resize(CTE::MAGcolumns,-1);
        for (unsigned int m = 0; m < unsigned_nMol; m++)
        {
            mol[m].MAG.clear();

            for (unsigned int j = 0; j < CTE::MAGmaxline; j++)  // j = 1 para considerar o cabecalho do arquivo
            {
                // Limpando vetor ligaçoes
                for (unsigned int idx = 0 ; idx < Ligacao.size() ; idx++) { Ligacao[idx] = -1; }

                //lendo a linha
                std::getline(file, line, '\n');
                std::stringstream iss(line);

                if (pulaCabecalho) pulaCabecalho = false;
                else
                {
                    for (unsigned int k = 0; k < CTE::MAGcolumns+2; k++)
                    {
                        //lendo a coluna
                        std::getline(iss, colum, ',');
                        
                        if (k > 1)
                        {
                            try 
                            {
                                Ligacao[k-2] = std::stoi(colum);
                            }
                            catch (const std::invalid_argument& e)
                            {
                                j = CTE::MAGmaxline;
                                k = CTE::MAGcolumns+2;
                            }
                            // if (m == 0)
                            // { 
                            //     std::cout << "\tk:" << k << " Ligacao:";
                            //     for (unsigned int idx = 0 ; idx < Ligacao.size() ; idx++) { std::cout << Ligacao[idx] << "," ; }
                            //     std::cout << std::endl;
                            // }
                        }
                    }
                    mol[m].MAG.push_back(Ligacao);
                }
            }
            Ligacao.resize(CTE::MAGcolumns,-1);
            mol[m].MAG.push_back(Ligacao);
        }
        // Fecha arquivo
        file.close();
    }
}

//Lê SOLex do arquivo
void Mistura::LeSOLex(std::string path, std::string filename)
{
    if (nMol < 1)std::cout << "Mistura::LeSOLex - Infome o numero de moleculas" << std::endl;
    else
    {
        path += "/" + filename;
        std::ifstream file(path);

        //Msg de Erro de abertura
        if (!file.is_open())
        {
            std::cout << "Mistura::LeSOLex - Nao foi possivel abrir " << path << std::endl;
            exit(EXIT_FAILURE);
        }

        //variaveis auxiliares
        std::string line, colum;
        for (unsigned int i = 0; i < mol.size(); i++)
        {
            mol[i].SOLex.clear();
            mol[i].SOLex.resize(CTE::Nucleosize); //dimensionando ao maximo de nucleos
            for (unsigned int j = 0; j < CTE::Nucleosize; j++) mol[i].SOLex[j].resize(CTE::SOLexpsize);    
        }

        //lendo as linhas
        unsigned int k = 0, i = 0;
        bool pulaCabecalho = true;
        while (std::getline(file, line))
        {
            
            if (pulaCabecalho) pulaCabecalho = false;
            else if (i < mol.size())
            {

                std::stringstream iss(line);

                for (unsigned int j = 0; j < CTE::SOLexpsize+2; j++)
                {
                    //lendo a coluna
                    std::getline(iss, colum, ',');
                   
                    if (j == 1)
                    {
                        try
                        {
                            std::string FimNucleo (1,colum.back());

                            // std::cout << "colum:" << colum << std::endl;
                            // std::cout << "colum.back():" << colum.back() << std::endl;
                            // std::cout << "FimNucleo:" << FimNucleo << std::endl;
                            // std::cout << "----------------------" << std::endl;

                            if (FimNucleo == "l") k = 0;
                            else k = std::stoll(FimNucleo);
                            
                        }
                        catch (const std::invalid_argument& e) 
                        {
                            // std::cout << "LeMSOLex: Entrou no catch" << std::endl;
                            k = 1;
                            j = CTE::SOLexpsize+2;
                            i++;
                        }
                    }
                    else if (j >= 2)
                    {
                        mol[i].SOLex[k][j-2] = std::stoll(colum);
                    }
                }
            }
        }
        file.close();
    }
}

//Lê composição do arquivo
void Mistura::LeComposicao(std::string path, std::string filename)
{
    //std::cout << "Le composicao" << std::endl << std::endl;
	path += "/" + filename;
	std::vector < std::vector < double > > M = LeMatriz(path,true); 
	double sumz = 0;
	double tol =  1e-6;

	z.resize(M.size(),1);
	for (unsigned int i = 0; i < M.size(); i++)
	{
		z(i,0) = M[i][1];
		sumz += z(i,0);
	}

	if (abs(1 - sumz) > tol)
	{
		for (unsigned int i = 0; i < M.size(); i++) z(i,0) = z(i,0)/sumz;
	}
}

//Imprime MAG no terminal
void Mistura::ImprimeMAG(int m)
{
    unsigned int inicio = 0;
    unsigned int fim = mol.size();
    if (m != -1)
    {// mostra apenas uma molécula 
        inicio = m;
        fim = m+1;
    }

    std::cout << "MAG:" << std::endl;
    std::cout << "mol,bond,Atom i,Atom j,Info,Core,Type" << std::endl;

    for (unsigned int m = inicio ; m < fim ; m++)
    {
        for (unsigned int i = 0 ; i < mol[m].MAG.size() ; i++)
        {//varre ligações
            std::cout << m << "," << i;
            for (unsigned int k = 0 ; k < mol[m].MAG[i].size(); k++)
            {//varre colunas
                std::cout << "," << mol[m].MAG[i][k];
            }
            if (mol[m].MAG[i][0] < 0) i = mol[m].MAG.size();
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

//Imprime SOLex no terminal
void Mistura::ImprimeSOLex(int m)
{
    unsigned int inicio = 0;
    unsigned int fim = nMol;
    if (m != -1)
    {// mostra apenas uma molécula 
        inicio = m;
        fim = m+1;
    }

    std::cout << "SOLex:" << std::endl;
    std::cout << "mol-core:A6,A4,A3,A2,N6,N5,N4a,N4n,N3a,N3m,N3n,N2n,N2a,N2m,N2f,N1a,N1m,N1n,Rp,Rm,Rn,Ra,br,br2,MEn,MEa,IHo,IHn,AAa,AAm,AAn,NS,RS,ANa,ANn,NN,RN,NO,RO,KO,Ni,VO,RCn,Arr" << std::endl;
    for (unsigned int i = inicio; i < fim; i++) 
    {
        for (int k = 0; k < mol[i].SOLex[0][CTE::SOLexpsize-1] + 1; k++)
        {
            std::cout << i << "-";
            if(k == 0) std::cout << "total" << ": ";
            else std::cout << "core" << k << ": ";
            for (unsigned int j = 0; j < mol[i].SOLex[k].size(); j++) std::cout << mol[i].SOLex[k][j] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

//Imprime MAGnode no terminal
void Mistura::ImprimeMAGnode(int m)
{
    unsigned int inicio = 0;
    unsigned int fim = nMol;
    if (m != -1)
    {// mostra apenas uma molécula 
        inicio = m;
        fim = m+1;
    }

    std::cout << "MAGnode:" << std::endl;
    for (unsigned int i = inicio; i < fim; i++) 
    {
        for (unsigned int j = 0; j < mol[i].MAGnode.size(); j++) 
        {
            std::cout << i << "-" << j << ": " << mol[i].MAGnode[j] << std::endl;
            if (mol[i].MAGnode[j] < 0) j = mol[i].MAGnode.size();
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//Imprime LigacoesPorAtomo no terminal
void Mistura::ImprimeLigacoesPorAtomo(int m)
{
    unsigned int inicio = 0;
    unsigned int fim = nMol;
    if (m != -1)
    {// mostra apenas uma molécula 
        inicio = m;
        fim = m+1;
    }

    std::cout << "LigacoesPorAtomo:" << std::endl;
    for (unsigned int i = inicio; i < fim; i++) 
    {
        for (unsigned int j = 0; j < mol[i].LigacoesPorAtomo.size(); j++) 
        {
            std::cout << i << "-" << j << ": ";
            for (unsigned int k = 0; k < mol[i].LigacoesPorAtomo[j].size(); k++) 
            {
                std::cout << mol[i].LigacoesPorAtomo[j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//Imprime a composição molar no temrinal
void Mistura::ImprimeComposicao(int m)
{
    unsigned int inicio = 0;
    unsigned int fim = z.rows();

    if (m != -1)
    {// mostra apenas uma molécula 
        inicio = m;
        fim = m+1;
    }

    std::cout << "Composicao:" << std::endl;
    for (unsigned int i = inicio; i < fim; i++) std::cout << i << ": " << z(i,0) << std::endl;
    std::cout << std::endl;
}

//Lê uma matriz 2D do arquivo "path" e retorna seu conteúdo num array 2D tipo double 
std::vector < std::vector < double > > Mistura::LeMatriz(std::string path, bool pulaCabecalho)
{
	std::vector < std::vector < double > > M; 
	std::ifstream file(path);

	//Msg de Erro de abertura
	if (!file.is_open())
	{
		std::cout << "Erro: nao foi possivel abrir " << path << std::endl;
		exit(EXIT_FAILURE);
	}

	//variaveis auxiliares
	std::string line, colum;
	unsigned int c = 0, l = 0;

	// encontrando a dimens�o da matriz de par�metros
	if (file.good())
	{
		//lendo as linhas
		while (std::getline(file, line))
		{
			// se estiver na primeira lilnha
			if (l == 0)
			{
				// Create a stringstream from line
				std::stringstream ss(line);

				//conta as colunas
				while (std::getline(ss, colum, ',')) c++;
			}

			//conta as linhas 
			l++;
		}
	}

	//voltando ao come�o do arrquivo
	file.clear();
	file.seekg(0);

    if (pulaCabecalho) // pulando cabeçalho quando necessário
    {
        l--;
        std::getline(file, line, '\n');
    }

	M.resize(l);
	for (unsigned int i = 0; i < l; i++)
	{
		M[i].resize(c);
		//lendo a linha
		std::getline(file, line, '\n');
		std::stringstream  iss(line);

		for (unsigned int j = 0; j < c; j++)
		{
			//lendo a coluna
			std::getline(iss, colum, ',');
			M[i][j] = std::stod(colum);
		}
	}
	// Fecha arquivo
	file.close();
	return M;
}

// Lê tipos de atomos das moleculas
void Mistura::CarregaMAGnodes(int molecula)
{   
    unsigned int inicio = 0;
    unsigned int fim = mol.size();

    if (molecula >= 0)
    {
        inicio = molecula;
        fim = molecula + 1;
    }

    for (unsigned int m = inicio ; m < fim ; m++)
    {
        ContaUltimoC(m);
        mol[m].MAGnode.resize(mol[m].MAG.size());

        int intMAGsize = mol[m].MAG.size();
        for (int i = 0 ; i < intMAGsize ; i++)
        {
            // varrendo átomos
            mol[m].MAGnode[i] = 0;
            for (unsigned int h = 0 ; h < mol[m].MAG.size() ; h++)
            {// varrrendo ligações
                if (i > mol[m].UltimoC) mol[m].MAGnode[i] = -1; //acabou a molécula
                else if ((mol[m].MAG[h][0] == i) || (mol[m].MAG[h][1] == i))
                {
                    if ((mol[m].MAG[h][2] >= 1000) && (mol[m].MAG[h][1] == i)) mol[m].MAGnode[i] = mol[m].MAG[h][2];    // encontramos um heteroátomo
                    else if  (mol[m].MAGnode[i] < 1000) mol[m].MAGnode[i]++;                                            // encontramos um carbono
                }
                if (mol[m].MAG[h+1][0] < 0) h = mol[m].MAG.size(); //Acabou a molécula
            }
        }
    }
}

//Verifica o ultimo carbono da molécula
void Mistura::ContaUltimoC(int molecula)
{
    unsigned int inicio = 0;
    unsigned int fim = mol.size();

    if (molecula >= 0)
    {
        inicio = molecula;
        fim = molecula + 1;
    }

    for (unsigned int m = inicio ; m < fim ; m++)
    {
        mol[m].UltimoC = 0;
        for (unsigned int i = 0 ; i < mol[m].MAG.size() ; i++)
        {
            if (mol[m].MAG[i][1] > mol[m].UltimoC) mol[m].UltimoC = mol[m].MAG[i][1];
            else if (mol[m].MAG[i][1] < 0) i = mol[m].MAG.size();  //Acabou a molécula
        }
    }
}

void Mistura::defSeed(int seed)
{
    Rnd.defSeed(seed);
}

int Mistura::BuscaLiga(int m, int Ci, int Cj)
{// Buscando a ligação que contem os Carbonos i e j

    int ligaProc = -1;
    for (unsigned int i = 0 ; i < mol[m].MAG.size() ; i++)
    {// varrendo todas as ligações da molécula j
        if (((mol[m].MAG[i][0] == Ci) && (mol[m].MAG[i][1] == Cj)) || ((mol[m].MAG[i][0] == Cj) && (mol[m].MAG[i][1] == Ci)))
        {
            ligaProc = i;
            i = mol[m].MAG.size(); //interrompe o for uma vez que encontrou a ligação
        }
    }
    return ligaProc;
}

// Carrega  ligações que contÊm o Atomo 'e' na molécula 'm'
void Mistura::CarregaLigacoesPorAtomo(int molecula, bool contaC)
{
    //std::cout << "DEBUG - CarregaLigacoesPorAtomo" << std::endl;
    unsigned int inicio = 0;
    unsigned int fim = mol.size();
    int hetero = 0;
    int conthetero = 0;

    if (molecula >= 0)
    {
        inicio = molecula;
        fim = molecula + 1;
    }
    for (unsigned int m = inicio ; m < fim ; m++)
    {
        if (contaC) ContaUltimoC(m);
        mol[m].LigacoesPorAtomo.clear();
        mol[m].LigacoesPorAtomo.resize(mol[m].UltimoC+2, std::vector<int>(CTE::NligasMaxPorAtomo,-1));
        mol[m].MAGnode.clear();
        mol[m].MAGnode.resize(mol[m].UltimoC+2,0);
        // if (m == 0) std::cout << m << " - UltimoC:" << mol[m].UltimoC << " " << mol[m].LigacoesPorAtomo.size() << "x" << mol[m].LigacoesPorAtomo[0].size() << ":" << mol[m].LigacoesPorAtomo[0][0] << "," << mol[m].LigacoesPorAtomo[0][1] << "," << mol[m].LigacoesPorAtomo[0][2] << "," << mol[m].LigacoesPorAtomo[0][3] << " mol[m].MAGnode[0]:" << mol[m].MAGnode[0] << std::endl; 
        for (unsigned int i = 0 ; i < mol[m].MAG.size() ; i++)
        {
            if (mol[m].MAG[i][0] < 0) i = mol[m].MAG.size();  //Acabou a molécula
            else
            {
                int Ai = mol[m].MAG[i][0];
                int Aj = mol[m].MAG[i][1];
                mol[m].LigacoesPorAtomo[Ai][mol[m].MAGnode[Ai]] = i;
                mol[m].LigacoesPorAtomo[Aj][mol[m].MAGnode[Aj]] = i;
                // if (m == 0) std::cout << "Ai:" << Ai << " Aj:" << Aj << " i:" << i << " mol[" << m << "].MAGnode[" << Ai << "]:" << mol[m].MAGnode[Ai] << " mol[" << m << "].MAGnode[" << Aj << "]:" << mol[m].MAGnode[Aj] << " mol[" << m << "].LigacoesPorAtomo[" << Ai << "][" << mol[m].MAGnode[Ai] << "]:" << mol[m].LigacoesPorAtomo[Ai][mol[m].MAGnode[Ai]] << std::endl;
                mol[m].MAGnode[Ai]++;
                mol[m].MAGnode[Aj]++;
            }
        }

        conthetero = 0;
        hetero = mol[m].SOLex[0][31] + mol[m].SOLex[0][32] + mol[m].SOLex[0][33] + mol[m].SOLex[0][34] + mol[m].SOLex[0][35] + mol[m].SOLex[0][36] + mol[m].SOLex[0][37] + mol[m].SOLex[0][38] + mol[m].SOLex[0][39]; 
        if (hetero > 0)
        {
            for (unsigned int i = 0 ; i < mol[m].MAG.size() ; i++)
            {
                if (mol[m].MAG[i][0] < 0 || conthetero >= hetero) i = mol[m].MAG.size();  //Acabou a molécula
                else
                {
                    if (mol[m].MAG[i][2] >= 1000) 
                    {
                        mol[m].MAGnode[mol[m].MAG[i][1]] = mol[m].MAG[i][2];
                        conthetero++;
                    }
                }
            }
        }
        mol[m].MAGnode[mol[m].MAGnode.size()-1] = -1;
    }
     //std::cout << "Termina - CarregaLigacoesPorAtomo" << std::endl;
}

void Mistura::CarregaSOLexEigen()
{
    SOLexEigen.resize(nMol,CTE::Mjoback_SOL.cols());
    for (unsigned int i = 0 ;  i < SOLexEigen.rows(); i++)
    {
        for (unsigned int j = 0 ;  j < SOLexEigen.cols(); j++)
        {
            SOLexEigen(i,j) = mol[i].SOLex[0][j];
        }
    }
}
