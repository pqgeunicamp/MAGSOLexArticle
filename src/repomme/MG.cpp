#include "MG.hpp"

namespace MG
{
    /////////////////////////////////////////////////////
    ////    FUNÇÕES PARA GRUPOS DE PRIMEIRA ORDEM    ////
    /////////////////////////////////////////////////////
    #pragma region Estruturas_1_ordem

    // Testa se o átomo 'A' é um carbono alifático de ponta de cadeia CH3 da molécula 'm' | -[CH3]
    bool CH3(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e1 = 0;
        bool e1tipo = 1;
        bool CH3cetona = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < std::vector < int > > GrupoCentona = GCH3CO(m, MX);

        // identificando se o CH2 está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][0] == A) CH3cetona = 1;
        }
        
        // O átomo 'A' tem que ser carbono primário 
        if (MX.mol[m].MAGnode[A] == 1 && CH3cetona == 0)
        {
            // buscando as ligações deste carbono 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // Pegando o primeiro elemento 'e1' que se liga a 'A'
            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e1 = MX.mol[m].MAG[Liga_A[0]][1];
            else e1 = MX.mol[m].MAG[Liga_A[0]][0];

            // testando se 'e1' é aromático
            e1tipo = pmul_e(m, e1, 2, MX);

            // o átomo 'e1' ligado a 'A' não pode ser aromático
            if (e1tipo == 0 && MX.mol[m].MAG[Liga_A[0]][2] == 999 && MX.mol[m].MAG[Liga_A[0]][4] == 0) _tipo = 1;
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono alifático não ramificado de meio de cadeia CH2 da molécula 'm' | -[CH2]-
    bool CH2(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e1 = 0;
        unsigned int e2 = 0;
        bool e1tipo = 1;
        bool e2tipo = 1;
        bool CH2cetona = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < std::vector < int > > GrupoCentona = GCH2COCH2(m, MX);


        // buscando as ligações deste carbono 'A'
        Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

        // identificando se o CH2 está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][0] == A) CH2cetona = 1;
        }
        
        if (MX.mol[m].MAGnode[A] == 2 && CH2cetona == 0 && ( ((MX.mol[m].MAG[Liga_A[0]][2] < 1000) || (MX.mol[m].MAG[Liga_A[0]][2] == 1002)) &&
        ((MX.mol[m].MAG[Liga_A[1]][2] < 1000) || (MX.mol[m].MAG[Liga_A[1]][2] == 1002)) )) 
        {  
            // As ligações deste carbono tem que ser saturadas alifáticas
            if ((MX.mol[m].MAG[Liga_A[0]][4] == 0) && (MX.mol[m].MAG[Liga_A[1]][4] == 0))
            {
                // Pegando o primeiro elemento que se liga a 'A' (elemento 'e1')
                if (MX.mol[m].MAG[Liga_A[0]][0] == A) e1 = MX.mol[m].MAG[Liga_A[0]][1];
                else e1 = MX.mol[m].MAG[Liga_A[0]][0];

                // testando se 'e1' é aromático
                e1tipo = pmul_e(m, e1, 2, MX);

                // o átomo 'e1' ligado a 'A' não pode ser aromático
                if (e1tipo == 0) 
                {
                    // Pegando o segundo elemento que se liga a 'A' (elemento 'e2')
                    if (MX.mol[m].MAG[Liga_A[1]][0] == A) e2 = MX.mol[m].MAG[Liga_A[1]][1];
                    else e2 = MX.mol[m].MAG[Liga_A[1]][0];

                    // testando se 'e2' é aromático
                    e2tipo = pmul_e(m, e2, 2, MX);

                    // o átomo 'e2' ligado a 'A' não pode ser aromático
                    if (e2tipo == 0) _tipo = 1;

                    //if (_tipo == 1 && m == 41) std::cout << "CH2 " << A << std::endl;
                }
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono alifático ramificado de meio de cadeia CH da molécula 'm' | >[CH]-
    bool CH(int m, int A, Mistura &MX)
    {
    // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e1 = 0;
        unsigned int e2 = 0;
        unsigned int e3 = 0;
        bool e1tipo = 1;
        bool e2tipo = 1;
        bool e3tipo = 1;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário 
        if (MX.mol[m].MAGnode[A] == 3 && CHpertenceKO(m, A, MX) == 0)
        {
            // buscando as ligações deste carbono 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // As ligações deste carbono tem que ser saturadas alifáticas
            if ((MX.mol[m].MAG[Liga_A[0]][4] == 0) && (MX.mol[m].MAG[Liga_A[1]][4] == 0) && (MX.mol[m].MAG[Liga_A[2]][4] == 0))
            {
                // Pegando o primeiro elemento que se liga a 'A' (elemento 'e1')
                if (MX.mol[m].MAG[Liga_A[0]][0] == A) e1 = MX.mol[m].MAG[Liga_A[0]][1];
                else e1 = MX.mol[m].MAG[Liga_A[0]][0];

                // testando se 'e1' é aromático
                e1tipo = pmul_e(m, e1, 2, MX);

                // o átomo 'e1' ligado a 'A' não pode ser aromático
                if (e1tipo == 0) 
                {
                    // Pegando o segundo elemento que se liga a 'A' (elemento 'e2')
                    if (MX.mol[m].MAG[Liga_A[1]][0] == A) e2 = MX.mol[m].MAG[Liga_A[1]][1];
                    else e2 = MX.mol[m].MAG[Liga_A[1]][0];

                    // testando se 'e2' é aromático
                    e2tipo = pmul_e(m, e2, 2, MX);

                    // o átomo 'e2' ligado a 'A' não pode ser aromático
                    if (e2tipo == 0) 
                    {
                        // Pegando o terceiro elemento que se liga a 'A' (elemento 'e3')
                        if (MX.mol[m].MAG[Liga_A[2]][0] == A) e3 = MX.mol[m].MAG[Liga_A[2]][1];
                        else e3 = MX.mol[m].MAG[Liga_A[2]][0];

                        // testando se 'e3' é aromático
                        e3tipo = pmul_e(m, e3, 2, MX);

                        // o átomo 'e3' ligado a 'A' não pode ser aromático
                        if (e3tipo == 0) _tipo = 1;
                    }
                }
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono aromático não substituído aCH da molécula 'm' | >[aCH]
    bool aCH (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // O átomo 'A' tem que ser carbono secundário 
        if (MX.mol[m].MAGnode[A] == 2) _tipo = pmul_e(m, A, 2, MX);
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono naftênico não substituído nCH2 da molécula 'm' | >[nCH2]
    bool nCH2 (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // O átomo 'A' tem que ser carbono secundário 
        if (MX.mol[m].MAGnode[A] == 2) _tipo = tal_e(m, A, 1, MX);
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono naftênico substituído nCH da molécula 'm'. Ligado a um alifático ou a um naftênico | >[nCH]-
    bool nCH (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // O átomo 'A' tem que ser carbono secundário 
        if (MX.mol[m].MAGnode[A] == 3 && CHpertenceKO(m, A, MX) == 0)
        {
            // testando se 'A' não é aromático e tem pelo menos 1 ligação naftênica
            if (pmul_e(m, A, 2, MX) == 0) _tipo = pmul_e(m, A, 1, MX);

            // if (_tipo == 1) std::cout << "nCH " << A << std::endl;
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono aromático fundido a outro anel aromático | >[aC]-aCHx-
    bool aCfaC (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int anel = 0;
        int e = 0;
        int anel_e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário 
        if (MX.mol[m].MAGnode[A] == 3)
        {
            // testando se não está lig à heteroátomo e se todas as ligações feitas por 'A' são aromáticas 
            if (toel_e(m, A, MX) == 1 && tal_e(m, A, 2, MX) == 1)
            {
                //identifica o anel de A
                anel = AnelPred(m, A, MX);

                // buscando as ligações deste átomo 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                //varrendo o vetor de ligações de A
                for (unsigned int i = 0; i < Liga_A.size(); i++)
                {
                    //identificando em qual posição está  o átomo ligado à A
                    if (MX.mol[m].MAG[Liga_A[i]][0]==A) e = MX.mol[m].MAG[Liga_A[i]][1];
                    else e = MX.mol[m].MAG[Liga_A[i]][0];

                    //identifica o anel de e
                    anel_e = AnelPred(m, e, MX);

                    //Identificando o carbono ligado à A que não pertence ao anel de A & se os dois anéis não são A6
                    if ((anel_e != anel) && (ContaLigaAnel(m, anel_e, MX) + ContaLigaAnel(m, anel, MX) < 12)) _tipo = 1;
    
                    //interrompe o "for" quando i é < 3, e a próxima lig =-1
                    if ((i < (Liga_A.size()-1)) && (Liga_A[i+1] == -1)) i = Liga_A.size();                          
                }
            }
            // if (_tipo == 1) std::cout << "aCfaC " << A << std::endl;
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono aromático fundido a um anel naftênico | >[aC]-nCHx-
    bool aCfnC (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário 
        if (MX.mol[m].MAGnode[A] == 3)
        {
            // testando se 'A' é aromático & não está ligado a heteroátomo & tem alguma ligação naftênica
            if ((pmul_e(m, A, 2, MX) == 1) && (toel_e(m, A, MX) == 1) && pmul_e(m, A, 1, MX) == 1)
            {
                // buscando as ligações deste átomo 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                //Varrendo o vetor de ligações do atomo A
                for (unsigned int i = 0; i < Liga_A.size();i++)
                {
                    //Testando se a ligação i é naftênica
                    if (MX.mol[m].MAG[Liga_A[i]][4] == 1)
                    {
                        //identificando em qual posição está  o átomo ligado à A
                        if (MX.mol[m].MAG[Liga_A[i]][0]==A) e = MX.mol[m].MAG[Liga_A[i]][1];
                        else e = MX.mol[m].MAG[Liga_A[i]][0];

                        //Identificar se não faz ligação aromática
                        if (pmul_e(m, e, 2, MX) == 0) _tipo = 1;
                    }
                    //interrompe o "for" quando i é < 3, e a próxima lig =-1
                    if ((i < (Liga_A.size()-1)) && (Liga_A[i+1] == -1)) i = Liga_A.size();
                }
            }
            //if (_tipo == 1) std::cout << "aCfnC " << A << std::endl;
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um carbono aromático ligado à outro carbono qq
    bool aCqq (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
    
        // O átomo 'A' tem que ser carbono terciário & aromático
        if ((MX.mol[m].MAGnode[A] == 3) && (pmul_e(m, A, 2, MX) == 1) && 
            (aCfaC(m, A, MX)+ aCfnC(m, A, MX) + aC_CH3(m, A, MX) + aC_CH2(m, A, MX) +
             aC_CH(m, A, MX) + aC_C(m, A, MX) + aCS(m, A, MX) + aCN(m, A, MX) + 
             aCO(m, A, MX) + aC_CO(m, A, MX) +  aldeido_aC_CO(m, A, MX) == 0))
        {
            _tipo = 1;
        }

        return _tipo;
    }

    // Testa se o átomo 'A' é um Enxofre ligado à aromático 
    bool aCS (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
    
        // O átomo 'A' tem que ser carbono terciário & aromático
        if ((MX.mol[m].MAGnode[A] == 3) && (pmul_e(m, A, 2, MX) == 1))
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            //varrendo o vetor de ligações de A
            for (unsigned int i = 0; i < Liga_A.size(); i++)
            { 
                if (MX.mol[m].MAG[Liga_A[i]][2] == 1000) _tipo = 1;

                //interrompe o "for" quando i é < 3, e a próxima lig =-1
                if ((i < (Liga_A.size()-1)) && (Liga_A[i+1] == -1)) i = Liga_A.size();       
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um Nitrogênio ligado à aromático. Ex carbazol arom - N - arom
    bool aCN (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);
    
        // O átomo 'A' tem que ser carbono terciário & aromático
        if ((MX.mol[m].MAGnode[A] == 3) && (pmul_e(m, A, 2, MX) == 1))
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            //varrendo o vetor de ligações de A
            for (unsigned int i = 0; i < Liga_A.size(); i++)
            { 
                if (MX.mol[m].MAG[Liga_A[i]][2] == 1001)
                {
                    e = MX.mol[m].MAG[Liga_A[i]][1]; //o próprio N
                    
                    // buscando as ligações deste átomo 'e', com o nitrogenio
                    Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

                    //identificando o átomo f (carbono), que tem que ser diferente do átomo A (carbono)
                    if (MX.mol[m].MAG[Liga_e[0]][0] == e) f = MX.mol[m].MAG[Liga_e[0]][1];
                    else f = MX.mol[m].MAG[Liga_e[0]][0];

                    if (f == A)
                    {
                        if (MX.mol[m].MAG[Liga_e[1]][0] == e) f = MX.mol[m].MAG[Liga_e[1]][1];
                        else f = MX.mol[m].MAG[Liga_e[1]][0];
                    }
                    
                    //identiicando que os anéis de A e f são distintos.
                    if (AnelPred(m, A, MX) != AnelPred(m, f, MX)) _tipo = 1;
                }
                //interrompe o "for" quando i é < 3, e a próxima lig =-1
                if ((i < (Liga_A.size()-1)) && (Liga_A[i+1] == -1)) i = Liga_A.size();       
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um Oxigênio ligado à aromático 
    bool aCO (int m, int A, Mistura &MX)
    {    
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
    
        // O átomo 'A' tem que ser carbono terciário & aromático
        if ((MX.mol[m].MAGnode[A] == 3) && (pmul_e(m, A, 2, MX) == 1))
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            //varrendo o vetor de ligações de A
            for (unsigned int i = 0; i < Liga_A.size(); i++)
            { 
                if (MX.mol[m].MAG[Liga_A[i]][2] == 1002) _tipo = 1;

                //interrompe o "for" quando i é < 3, e a próxima lig =-1
                if ((i < (Liga_A.size()-1)) && (Liga_A[i+1] == -1)) i = Liga_A.size();       
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um carbono aromático está ligado a um metil | [>aC-CH3]
    bool aC_CH3(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário
        if (MX.mol[m].MAGnode[A] == 3)
        {
            // testando se 'A' é aromático
            if (pmul_e(m, A, 2, MX) == 1)
            {
                // buscando as ligações deste carbono 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                // varrendo ligações de 'A'
                for (unsigned int L = 0 ; L < Liga_A.size() ; L++)
                {
                    // Pegando o elemento 'e' que se liga a 'A'
                    if (MX.mol[m].MAG[Liga_A[L]][0] == A) e = MX.mol[m].MAG[Liga_A[L]][1];
                    else e = MX.mol[m].MAG[Liga_A[L]][0];

                    // Se o Carbono 'e' for primário e saturado conta 
                    if (MX.mol[m].MAGnode[e] == 1) _tipo = 1;

                    // se não existe a próxima ligação termina o 'for'
                    if (Liga_A[L+1] < 0) L = Liga_A.size();
                }
            } 
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um carbono aromático está ligado a um carbono alifático secundário | [>aC-CH2-]
    bool aC_CH2(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
    
        // O átomo 'A' tem que ser carbono terciário
        if (MX.mol[m].MAGnode[A] == 3 && aCpertence_aCCH2aC(m, A, MX) == 0)
        {
            // testando se 'A' é aromático
            if (pmul_e(m, A, 2, MX) == 1)
            {
                // buscando as ligações deste carbono 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                // varrendo ligações de 'A'
                for (unsigned int L = 0 ; L < Liga_A.size() ; L++)
                {
                    // Pegando o elemento 'e' que se liga a 'A'
                    if (MX.mol[m].MAG[Liga_A[L]][0] == A) e = MX.mol[m].MAG[Liga_A[L]][1];
                    else e = MX.mol[m].MAG[Liga_A[L]][0];
                    
                    // Se o Carbono 'e' for secundario, não forma grupo centona (CH2-C=O) e é saturado 
                    if (MX.mol[m].MAGnode[e] == 2 && CH2pertenceKO(m, e, MX) == 0) _tipo = tal_e(m, e, 0, MX);
                    
                    // se não existe a próxima ligação termina o 'for'
                    if (Liga_A[L+1] < 0) L = Liga_A.size();
                }
            } 
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um carbono aromático está ligado a um carbono alifático terciário | [>aC-CH<]
    bool aC_CH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário
        if (MX.mol[m].MAGnode[A] == 3)
        {
            // testando se 'A' é aromático
            if (pmul_e(m, A, 2, MX) == 1)
            {
                // buscando as ligações deste carbono 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                // varrendo ligações de 'A'
                for (unsigned int L = 0 ; L < Liga_A.size() ; L++)
                {
                    // Pegando o elemento 'e' que se liga a 'A'
                    if (MX.mol[m].MAG[Liga_A[L]][0] == A) e = MX.mol[m].MAG[Liga_A[L]][1];
                    else e = MX.mol[m].MAG[Liga_A[L]][0];

                    // Se o Carbono 'e' for Terciário e saturado 
                    if (MX.mol[m].MAGnode[e] == 3 && CHpertenceKO(m, e, MX) == 0) _tipo = tal_e(m, e, 0, MX);

                    // se não existe a próxima ligação termina o 'for'
                    if (Liga_A[L+1] < 0) L = Liga_A.size();
                }
            } 
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um carbono aromático está ligado a um carbono alifático quaternário | [>aC-C<-]
    bool aC_C(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono terciário
        if (MX.mol[m].MAGnode[A] == 3)
        {
            // testando se 'A' é aromático
            if (pmul_e(m, A, 2, MX) == 1)
            {
                // buscando as ligações deste carbono 'A'
                Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

                // varrendo ligações de 'A'
                for (unsigned int L = 0 ; L < Liga_A.size() ; L++)
                {
                    // Pegando o elemento 'e' que se liga a 'A'
                    if (MX.mol[m].MAG[Liga_A[L]][0] == A) e = MX.mol[m].MAG[Liga_A[L]][1];
                    else e = MX.mol[m].MAG[Liga_A[L]][0];

                    // Se o Carbono 'e' for quaternário e saturado 
                    if (MX.mol[m].MAGnode[e] == 4) _tipo = tal_e(m, e, 0, MX);

                    // se não existe a próxima ligação termina o 'for'
                    if (Liga_A[L+1] < 0) L = Liga_A.size();
                }
            } 
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um carbono aromático está ligado a um carbono da cetona [>aC-C=O<-]
    bool aC_CO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GaC_CO(m, MX);   

        // identificando se o aC está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][0] == A) _tipo = 1;
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um Nitrogênio em anel aromático aN (Piridina) da molécula 'm' | =aC-[aN]=aC-
    bool aN (int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // O átomo 'A' tem que ser Nitrogêneo 
        if (MX.mol[m].MAGnode[A] == 1001) _tipo = tal_e(m, A, 2, MX);
        return _tipo;
    }

    // Testa se o átomo 'A' é um CH2SH ligado a cadeia
    bool CH2SH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        int liga_ef = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser Enxofre 
        if (MX.mol[m].MAGnode[A] == 1000)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // identificando a posição do átomo "e ligado a "A"
            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

            // testando se "e" é um carbono secundário & saturado
            if (MX.mol[m].MAGnode[e] == 2 && tal_e(m, e, 0, MX) ==1)
            {
                // buscando as ligações deste átomo 'e'
                Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

                // Identificando a outra ligação que "e" faz sem ser a com o "A", ie ligação com "f"
                if (Liga_e[0] == Liga_A[0]) liga_ef = Liga_e[1];
                else liga_ef = Liga_e[0];
            
                if (MX.mol[m].MAG[liga_ef][0] == e) f = MX.mol[m].MAG[liga_ef][1];
                else f = MX.mol[m].MAG[liga_ef][0];

                // Identificando se o "f" não é aromático
                if (pmul_e(m, f, 2, MX) == 0) _tipo = 1;
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um SH ligado a um grupamento prioritário
    bool SH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        int liga_ef = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser Enxofre 
        if (MX.mol[m].MAGnode[A] == 1000)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // identificando a posição do átomo "e" ligado a "A"
            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

            // testando se "e" é um carbono secundário & saturado
            if (MX.mol[m].MAGnode[e] == 2 && tal_e(m, e, 0, MX) ==1)
            {
                // buscando as ligações deste átomo 'e'
                Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

                // Identificando a outra ligação que "e" faz sem ser a com o "A", ie ligação com "f"
                if (Liga_e[0] == Liga_A[0]) liga_ef = Liga_e[1];
                else liga_ef = Liga_e[0];
            
                if (MX.mol[m].MAG[liga_ef][0] == e) f = MX.mol[m].MAG[liga_ef][1];
                else f = MX.mol[m].MAG[liga_ef][0];

                // Identificando se o "f" é aromático
                if (pmul_e(m, f, 2, MX) == 1) _tipo = 1;
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um NH2 ligado a um grupamento prioritário
    bool NH2(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        int liga_ef = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser Nitrogênio 
        if (MX.mol[m].MAGnode[A] == 1001)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // identificando a posição do átomo "e" ligado a "A"
            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

            // testando se "e" é um carbono secundário & saturado
            if (MX.mol[m].MAGnode[e] == 2 && tal_e(m, e, 0, MX) ==1)
            {
                // buscando as ligações deste átomo "e"
                Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

                // Identificando a outra ligação que "e" faz sem ser a com o "A", ie ligação identificada como "f"
                if (Liga_e[0] == Liga_A[0]) liga_ef = Liga_e[1];
                else liga_ef = Liga_e[0];
            
                if (MX.mol[m].MAG[liga_ef][0] == e) f = MX.mol[m].MAG[liga_ef][1];
                else f = MX.mol[m].MAG[liga_ef][0];

                // Identificando se o "f" é aromático
                if (pmul_e(m, f, 2, MX) == 1) _tipo = 1;
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um CH2NH2 ligado a cadeia
    bool CH2NH2(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        int liga_ef = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser Nitrogênio 
        if (MX.mol[m].MAGnode[A] == 1001)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);
            if (MX.mol[m].MAG[Liga_A[0]][0]==A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

        
            // testando se "e" é um carbono secundário & saturado
            if (MX.mol[m].MAGnode[e] == 2 && tal_e(m, e, 0, MX) ==1)
            {
                // buscando as ligações deste átomo 'e'
                Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

                // Identificando a outra ligação que "e" faz sem ser a com o "A", ie ligação com "f"
                if (Liga_e[0] == Liga_A[0]) liga_ef = Liga_e[1];
                else liga_ef = Liga_e[0];
            
                if (MX.mol[m].MAG[liga_ef][0] == e) f = MX.mol[m].MAG[liga_ef][1];
                else f = MX.mol[m].MAG[liga_ef][0];

                // Identificando se o "f" não é aromático
                if (pmul_e(m, f, 2, MX) == 0) _tipo = 1;
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' é um OH ligado a cadeia
    bool OH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        if (MX.mol[m].MAGnode[A] == 1002)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

            //o átomo ligado ao "O" não pode ser aromático e o carbono deve ser secundário (o que exclui formação de cetona)
            if (pmul_e(m, e, 2, MX) == 0) _tipo = 1;    
        }
        return _tipo;
    }

    // Grupo cetona ligado a CH3 
    bool CH3CO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GCH3CO(m, MX);   

        // identificando se o CH3 está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = GrupoCentona.size();
        }
        return _tipo;
    }

    // Grupo cetona ligado a CH2 
    bool CH2CO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GCH2COCH2(m, MX);   

        // identificando se o CH2 está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = GrupoCentona.size();
        }
        return _tipo;
    }

    // Grupo cetona ligado a CH 
    bool CHCO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GCHCO(m, MX);   

        // identificando se o C  centônico (C=O) é desse tipo de cetona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = GrupoCentona.size();
        }
        return _tipo;
    }

    // Grupo R-C=O da função aldeído
    bool aldeido_CO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        int e = 0;
        int f = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        //encontrando o oxigenio do grupo C=O do aldeído
        if (MX.mol[m].MAGnode[A] == 1003)
        {
            // buscando as ligações deste átomo 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            if (MX.mol[m].MAG[Liga_A[0]][0] == A) e = MX.mol[m].MAG[Liga_A[0]][1];
            else e = MX.mol[m].MAG[Liga_A[0]][0];

            //o átomo de C ligado ao "O" não pode ser aromático e o carbono deve ser secundário para ser aldeído
            if (MX.mol[m].MAGnode[e] == 2 && pmul_e(m, e, 2, MX) == 0) 
            {
                // buscando as ligações deste átomo 'e' que é o carbono ligado ao O
                Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);
                
                // encontrando o outro carbono ligado ao  "e"           
                if (MX.mol[m].MAG[Liga_e[0]][0] == e) f = MX.mol[m].MAG[Liga_e[0]][1];
                else f = MX.mol[m].MAG[Liga_e[0]][0];
                
                //O carbono ligado ao "e" não pode ser aromático
                if (pmul_e(m, f, 2, MX) == 0) _tipo = 1;
            }
        }
        return _tipo;
    }    

    // Identifica se o aC é do Grupo aC-C=O da função aldeído
    bool aldeido_aC_CO(int m, int A, Mistura &MX)
    {
        
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > Grupo = Galdeido_aC_CO (m, MX);   

        // identificando se o aC está ligado a um grupo aldeído
        for (unsigned int i = 0; i < Grupo.size(); i++)
        {
            if (Grupo[i][0] == A) _tipo = 1;
            if (_tipo == 1) i = Grupo.size();
        }
        return _tipo;
    }    

    // Testa se o átomo 'A' é um carbono quaternário alifático ramificado de meio de cadeia C da molécula 'm' | >[C]<
    bool C4rio(int m, int A, Mistura &MX)
    {
    // variáveis auxiliares da função
        bool _tipo = 0;
        unsigned int e1 = 0;
        unsigned int e2 = 0;
        unsigned int e3 = 0;
        unsigned int e4 = 0;
        bool e1tipo = 1;
        bool e2tipo = 1;
        bool e3tipo = 1;
        bool e4tipo = 1;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);

        // O átomo 'A' tem que ser carbono quaternário
        if (MX.mol[m].MAGnode[A] == 4)
        {
            // buscando as ligações deste carbono 'A'
            Liga_A = MX.mol[m].LigacoesPorAtomo[A];//BuscaLigas_e(m, A, MX);

            // As ligações deste carbono tem que ser saturadas alifáticas
            if ((MX.mol[m].MAG[Liga_A[0]][4] == 0) && (MX.mol[m].MAG[Liga_A[1]][4] == 0) && (MX.mol[m].MAG[Liga_A[2]][4] == 0) && (MX.mol[m].MAG[Liga_A[3]][4] == 0))
            {
                // Pegando o primeiro elemento que se liga a 'A' (elemento 'e1')
                if (MX.mol[m].MAG[Liga_A[0]][0] == A) e1 = MX.mol[m].MAG[Liga_A[0]][1];
                else e1 = MX.mol[m].MAG[Liga_A[0]][0];

                // testando se 'e1' é aromático
                e1tipo = pmul_e(m, e1, 2, MX);

                // o átomo 'e1' ligado a 'A' não pode ser aromático
                if (e1tipo == 0)
                {
                    // Pegando o segundo elemento que se liga a 'A' (elemento 'e2')
                    if (MX.mol[m].MAG[Liga_A[1]][0] == A) e2 = MX.mol[m].MAG[Liga_A[1]][1];
                    else e2 = MX.mol[m].MAG[Liga_A[1]][0];

                    // testando se 'e2' é aromático
                    e2tipo = pmul_e(m, e2, 2, MX);

                    // o átomo 'e2' ligado a 'A' não pode ser aromático
                    if (e2tipo == 0)
                    {
                        // Pegando o terceiro elemento que se liga a 'A' (elemento 'e3')
                        if (MX.mol[m].MAG[Liga_A[2]][0] == A) e3 = MX.mol[m].MAG[Liga_A[2]][1];
                        else e3 = MX.mol[m].MAG[Liga_A[2]][0];

                        // testando se 'e3' é aromático
                        e3tipo = pmul_e(m, e3, 2, MX);

                        // o átomo 'e2' ligado a 'A' não pode ser aromático
                        if (e3tipo == 0)
                        {
                            // Pegando o quarto elemento que se liga a 'A' (elemento 'e4')
                            if (MX.mol[m].MAG[Liga_A[3]][0] == A) e4 = MX.mol[m].MAG[Liga_A[3]][1];
                            else e4 = MX.mol[m].MAG[Liga_A[3]][0];

                            // testando se 'e4' é aromático
                            e4tipo = pmul_e(m, e4, 2, MX);

                            // o átomo 'e4' ligado a 'A' não pode ser aromático
                            if (e4tipo == 0) _tipo = 1;
                        }
                    }
                }
            }
        }
        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um CH2 está ligado com ligacao dupla a um CH | [CH2=CH-]
    bool CH2dCH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        int nrLigInsaturadas = nrlt_e(m, A, 3, MX);

        if (nrLigInsaturadas == 1)
        {
            // identificando se o CH está ligado a um grupo CH2=CH-R
            for (unsigned int i = 0; i < MX.mol[m].mGCH2dCH.size(); i++)
            {
                if (MX.mol[m].mGCH2dCH[i][1] == A ) _tipo = 1;
                if (_tipo == 1) i = MX.mol[m].mGCH2dCH.size();
            }
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um CH2 está ligado com ligacao dupla a um CH | [CH2=C<]
    bool CH2dC(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        int nrLigInsaturadas = nrlt_e(m, A, 3, MX);

        if (nrLigInsaturadas == 1)
        {
            // identificando se o C está ligado a um grupo CH2=C<
            for (unsigned int i = 0; i < MX.mol[m].mGCH2dC.size(); i++)
            {
                if (MX.mol[m].mGCH2dC[i][1] == A) _tipo = 1;
                if (_tipo == 1) i = MX.mol[m].mGCH2dC.size();
            }
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um CH está ligado com ligacao dupla a um CH | [-CH=CH-]
    bool CHdCH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        bool CH2dCdCH = 0;
        bool CH2dCdC = 0;
        bool CHdCdCH = 0;
        
        //identificando se o atomo A é o CH2 do grupo (CH2=C=CH)-R
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdCH.size(); i++)
        {
            if ( MX.mol[m].mGCH2dCdCH[i][0] == A || MX.mol[m].mGCH2dCdCH[i][1] == A || MX.mol[m].mGCH2dCdCH[i][2] == A) CH2dCdCH = 1;
        }

        // identificando se o atomo A é do grupo (CH2=C=C<)
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdC.size(); i++)
        {
            if ( MX.mol[m].mGCH2dCdC[i][0] == A || MX.mol[m].mGCH2dCdC[i][1] == A || MX.mol[m].mGCH2dCdC[i][2] == A) CH2dCdC = 1;
        }

        // identificando se o atomo A é do grupo R-(CH=C=CH)-R
        for (unsigned int i = 0; i < MX.mol[m].mGCHdCdCH.size(); i++)
        {
            if ( MX.mol[m].mGCHdCdCH[i][0] == A || MX.mol[m].mGCHdCdCH[i][1] == A || MX.mol[m].mGCHdCdCH[i][2] == A) CHdCdCH = 1;
        }

        if (CH2dCdCH == 0 && CHdCdCH == 0 && CH2dCdC == 0)
        {
            // identificando se o primeiro CH está ligado a um grupo R-CH=CH-R
            for (unsigned int i = 0; i < MX.mol[m].mGCHdCH.size(); i++)
            {
                if (MX.mol[m].mGCHdCH[i][0] == A) _tipo = 1;
                if (_tipo == 1) i = MX.mol[m].mGCHdCH.size();
            }
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um CH2 está ligado com ligacao dupla a um CH | [CH2=C<]
    bool CHdC(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        int nrLigInsaturadas = nrlt_e(m, A, 3, MX);

        if (nrLigInsaturadas == 1)
        {
            // identificando se o C está ligado a um grupo CH2=C<
            for (unsigned int i = 0; i < MX.mol[m].mGCHdC.size(); i++)
            {
                if (MX.mol[m].mGCHdC[i][1] == A) _tipo = 1;
                if (_tipo == 1) i = MX.mol[m].mGCHdC.size();
            }
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe uma dupla de átomos onde um CH está ligado com ligacao dupla a um CH | [-CH=CH-]
    bool CdC(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        bool CH2dCdCH = 0;
        bool CH2dCdC = 0;
        bool CHdCdCH = 0;
        
        //identificando se o atomo A é o CH2 do grupo (CH2=C=CH)-R
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdCH.size(); i++)
        {
            if ( MX.mol[m].mGCH2dCdCH[i][0] == A || MX.mol[m].mGCH2dCdCH[i][1] == A || MX.mol[m].mGCH2dCdCH[i][2] == A) CH2dCdCH = 1;
        }

        // identificando se o atomo A é do grupo (CH2=C=C<)
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdC.size(); i++)
        {
            if ( MX.mol[m].mGCH2dCdC[i][0] == A || MX.mol[m].mGCH2dCdC[i][1] == A || MX.mol[m].mGCH2dCdC[i][2] == A) CH2dCdC = 1;
        }

        // identificando se o atomo A é do grupo R-(CH=C=CH)-R
        for (unsigned int i = 0; i < MX.mol[m].mGCHdCdCH.size(); i++)
        {
            if ( MX.mol[m].mGCHdCdCH[i][0] == A || MX.mol[m].mGCHdCdCH[i][1] == A || MX.mol[m].mGCHdCdCH[i][2] == A) CHdCdCH = 1;
        }

        if (CH2dCdCH == 0 && CHdCdCH == 0 && CH2dCdC == 0)
        {
            // identificando se o primeiro CH está ligado a um grupo R-CH=CH-R
            for (unsigned int i = 0; i < MX.mol[m].mGCdC.size(); i++)
            {
                if (MX.mol[m].mGCdC[i][0] == A) _tipo = 1;
                if (_tipo == 1) i = MX.mol[m].mGCdC.size();
            }
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe um grupo de átomos do tipo  [CH2=C=CH-]
    bool CH2dCdCH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // identificando se o =C= está ligado a um grupo [CH2=C=CH-]
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdCH.size(); i++)
        {
            if (MX.mol[m].mGCH2dCdCH[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = MX.mol[m].mGCH2dCdCH.size();
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe um grupo de átomos do tipo  [CH2=C=C<]
    bool CH2dCdC(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // identificando se o =C= está ligado a um grupo [CH2=C=C<]
        for (unsigned int i = 0; i < MX.mol[m].mGCH2dCdC.size(); i++)
        {
            if (MX.mol[m].mGCH2dCdC[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = MX.mol[m].mGCH2dCdC.size();
        }

        return _tipo;
    }

    // Testa se o átomo 'A' compõe um grupo de átomos do tipo  [-CH=C=CH-]
    bool CHdCdCH(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;

        // identificando se o =C= está ligado a um grupo [CH2=C=CH-]
        for (unsigned int i = 0; i < MX.mol[m].mGCHdCdCH.size(); i++)
        {
            if (MX.mol[m].mGCHdCdCH[i][1] == A) _tipo = 1;
            if (_tipo == 1) i = MX.mol[m].mGCHdCdCH.size();
        }

        return _tipo;
    }

    // Identificação a estrutura do grupo (CH2=CH)-R
    std::vector < std::vector < int > >  GCH2dCH(int m, Mistura &MX)
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(2);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            //encontrando a ligacao dupla, Condicoes:
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=CH)-R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 2) ou (Ci == 2 e Cj == 1) 
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                bool c1 = false; 
                bool c2 = false;
                int CH2, CH;
                
                //if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 1 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2 )) 
                if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 1 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2 )) 
                { 
                    c1 = true;
                    CH2 = MX.mol[m].MAG[j][0]; //Ci
                    CH = MX.mol[m].MAG[j][1]; //Cj
                }
                else if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 1 )) 
                {
                    c2 = true;
                    CH2 = MX.mol[m].MAG[j][1]; //Cj
                    CH = MX.mol[m].MAG[j][0]; //Ci
                }

                if (c1 == true || c2 == true)
                {
                    //"if" dentro de outro "if"
                    if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1)
                    {
                        resposta[0][0] = CH2; 
                        resposta[0][1] = CH; 
                    } else {
                        std::vector < int > _newMG = { CH2, CH };
                        // push_back em vetor bidimensional _newMG na matriz "resposta". 
                        // std::move e mais eficiente computacionalmente e usa menos linhas.
                        resposta.push_back(std::move(_newMG));
                    }
                }
            }
        }

        return resposta;
    }

    // Identificação a estrutura do grupo (CH2=C<)
    std::vector < std::vector < int > >  GCH2dC(int m, Mistura &MX)
    {
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(2);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            bool c1 = false; 
            bool c2 = false;
            int CH2, C;

            //encontrando a ligacao dupla, condicoes
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=C<) entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 3) ou (Ci == 3 e Cj == 1)  
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 1 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 3 )) 
                { 
                    c1 = true;
                    CH2 = MX.mol[m].MAG[j][0];
                    C = MX.mol[m].MAG[j][1];
                }
                else if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 3 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 1 )) 
                {
                    c2 = true;
                    CH2 = MX.mol[m].MAG[j][1];
                    C = MX.mol[m].MAG[j][0];
                }

                if (c1 == true || c2 == true)
                {
                    //"if" dentro de outro "if"
                    if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1)
                    {
                        resposta[0][0] = CH2;
                        resposta[0][1] = C; 
                    } else {
                        std::vector < int > _newMG = { CH2, C };
                        // push_back em vetor bidimensional _newMG na matriz "resposta". 
                        // std::move e mais eficiente computacionalmente e usa menos linhas.
                        resposta.push_back(std::move(_newMG));
                    }
                }
            }
        }

        return resposta;
    }

    // Identificação a estrutura do grupo (CH=C<)
    std::vector < std::vector < int > >  GCHdC(int m, Mistura &MX)
    {
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(2);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            bool c1 = false; 
            bool c2 = false;
            int CH, C;

            //encontrando a ligacao dupla, condicoes
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=C<) entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 3) ou (Ci == 3 e Cj == 1)  
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 3 )) 
                { 
                    c1 = true;
                    CH = MX.mol[m].MAG[j][0];
                    C = MX.mol[m].MAG[j][1];
                }
                else if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 3 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2 )) 
                {
                    c2 = true;
                    CH = MX.mol[m].MAG[j][1];
                    C = MX.mol[m].MAG[j][0];
                }

                if (c1 == true || c2 == true)
                {
                    //"if" dentro de outro "if"
                    if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1)
                    {
                        resposta[0][0] = CH;
                        resposta[0][1] = C; 
                    } else {
                        std::vector < int > _newMG = { CH, C };
                        // push_back em vetor bidimensional _newMG na matriz "resposta". 
                        // std::move e mais eficiente computacionalmente e usa menos linhas.
                        resposta.push_back(std::move(_newMG));
                    }
                }
            }
        }

        return resposta;
    }

    // Identificação a estrutura do grupo R-(CH=CH)-R
    std::vector < std::vector < int > >  GCHdCH(int m, Mistura &MX)
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(2);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            //encontrando a ligacao dupla, condicoes
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que R-(CH=CH)-R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == Cj == 2)  
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999 && ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2 )))
            {
                if (resposta.size() == unsigned(1) && resposta[0][0] == -1  && resposta[0][1] == -1 )
                {
                    resposta[0][0] = MX.mol[m].MAG[j][0];
                    resposta[0][1] = MX.mol[m].MAG[j][1];
                } else {
                    std::vector < int > _newMG = {MX.mol[m].MAG[j][0], MX.mol[m].MAG[j][1]};
                    // push_back em vetor bidimensional _newMG na matriz "resposta". 
                    // std::move e mais eficiente computacionalmente e usa menos linhas.
                    resposta.push_back(std::move(_newMG)); 
                }
            } 
        }

        return resposta;
    }

    // Identificação a estrutura do grupo R>(C=C)<R
    std::vector < std::vector < int > > GCdC(int m, Mistura &MX)
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(2);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            //encontrando a ligacao dupla, condicoes
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que R>(C=C)<R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == Cj == 2)  
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999 && ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 3 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 3 )))
            {
                if (resposta.size() == unsigned(1) && resposta[0][0] == -1  && resposta[0][1] == -1 )
                {
                    resposta[0][0] = MX.mol[m].MAG[j][0];
                    resposta[0][1] = MX.mol[m].MAG[j][1];
                } else {
                    std::vector < int > _newMG = {MX.mol[m].MAG[j][0], MX.mol[m].MAG[j][1]};
                    // push_back em vetor bidimensional _newMG na matriz "resposta". 
                    // std::move e mais eficiente computacionalmente e usa menos linhas.
                    resposta.push_back(std::move(_newMG)); 
                }
            } 
        }

        return resposta;
    }

    std::vector < std::vector < int > >  GCH2dCdCH(int m, Mistura &MX) 
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(3);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        resposta[0][2] = -1;

        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {
            int _CH2 = -1; 
            int _C = -1;
            int _CH = -1;

            //encontrando a PRIMEIRA ligacao dupla (CH2=C), Condicoes:
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=CH)-R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 2) ou (Ci == 2 e Cj == 1)
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] + MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 3 )
                {
                    if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 1) {
                        _CH2 = MX.mol[m].MAG[j][0];
                    } else if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2) {
                        _C = MX.mol[m].MAG[j][0];
                    }

                    if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 1) {
                        _CH2 = MX.mol[m].MAG[j][1];
                    } else if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2) {
                        _C = MX.mol[m].MAG[j][1];
                    }
                    
                    // buscando as ligações deste átomo '=C='
                    //Liga_C = BuscaLigas_e(m, _C);
                    Liga_C = MX.mol[m].LigacoesPorAtomo[_C];
                    
                    for (unsigned int li = 0; li < Liga_C.size(); li++)
                    { 
                        if ( (Liga_C[li] > -1) && (Liga_C[li] != signed(j)) && (MX.mol[m].MAG[Liga_C[li]][2] == 999) && (MX.mol[m].MAG[Liga_C[li]][4] == 3)) 
                        {
                            //verifica se e o Carbono (=C=) em ambas ligacoes e o mesmo, se for escolhe o seu companheiro da ligacao (CH).
                            if ((MX.mol[m].MAG[Liga_C[li]][0] == _C) && (MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]][1]] == 2)) {
                                // MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]]
                                _CH = MX.mol[m].MAG[Liga_C[li]][1]; 
                            } else if ((MX.mol[m].MAG[Liga_C[li]][1] == _C) && (MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]][0]] == 2)) {
                                _CH = MX.mol[m].MAG[Liga_C[li]][0];
                            }
                            
                            if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1 && resposta[0][2] == -1 && _CH != -1)
                            {
                                resposta[0][0] = _CH2;
                                resposta[0][1] = _C; 
                                resposta[0][2] = _CH;
                            } else if (_CH != -1) {
                                std::vector < int > newMG = { _CH2, _C, _CH }; 
                                resposta.push_back(std::move(newMG));
                            }
                        }
                    }  
                } 
            }
        }

        return resposta;
    }

    std::vector < std::vector < int > >  GCH2dCdC(int m, Mistura &MX) 
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(3);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        resposta[0][2] = -1;

        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            int _CH2 = -1; 
            int _C1 = -1;
            int _C2 = -1;

            //encontrando a PRIMEIRA ligacao dupla (CH2=C), Condicoes:
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=CH)-R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 2) ou (Ci == 2 e Cj == 1)
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] + MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 3 )
                {
                    if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 1) {
                        _CH2 = MX.mol[m].MAG[j][0];
                    } else if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2) {
                        _C1 = MX.mol[m].MAG[j][0];
                    }

                    if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 1) {
                        _CH2 = MX.mol[m].MAG[j][1];
                    } else if (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2) {
                        _C1 = MX.mol[m].MAG[j][1];
                    }
                    
                    // buscando as ligações deste átomo '=C='
                    //Liga_C = BuscaLigas_e(m, _C1);
                    Liga_C = MX.mol[m].LigacoesPorAtomo[_C1];
                    
                    for (unsigned int li = 0; li < Liga_C.size(); li++)
                    { 
                        if ( (Liga_C[li] > -1) && (Liga_C[li] != signed(j)) && (MX.mol[m].MAG[Liga_C[li]][2] == 999) && (MX.mol[m].MAG[Liga_C[li]][4] == 3)) 
                        {
                            //verifica se e o Carbono (=C1=) em ambas ligacoes e o mesmo, se for escolhe o seu companheiro da ligacao (C2).
                            if ((MX.mol[m].MAG[Liga_C[li]][0] == _C1) && (MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]][1]] == 3)) {
                                _C2 = MX.mol[m].MAG[Liga_C[li]][1]; 
                            } else if ((MX.mol[m].MAG[Liga_C[li]][1] == _C1) && (MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]][0]] == 3)) {
                                _C2 = MX.mol[m].MAG[Liga_C[li]][0];
                            }
                            
                            if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1 && resposta[0][2] == -1 && _C2 != -1 && _C1 != -1)
                            {
                                resposta[0][0] = _CH2;
                                resposta[0][1] = _C1; 
                                resposta[0][2] = _C2;
                            } else if (_C2 != -1 && _C1 != -1) {
                                std::vector < int > newMG = { _CH2, _C1, _C2 }; 
                                resposta.push_back(std::move(newMG));
                            }
                        }
                    }  
                } 
            }
        }

        return resposta;
    }

    std::vector < std::vector < int > >  GCHdCdCH(int m, Mistura &MX)
    {

        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        resposta.resize(1);
        resposta[0].resize(3);
        resposta[0][0] = -1;
        resposta[0][1] = -1;
        resposta[0][2] = -1;

        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++)
        {
            int CH_A, C, CH_B;

            //encontrando a PRIMEIRA ligacao dupla (CH2=C), Condicoes:
            //1. MX.mol[m].MAG[j][4] == 3   ---> Tipo de ligacao == Insaturada alifática
            //2. MX.mol[m].MAG[j][2] == 999  ---> Informação de posição da ligação == ligação **entre carbonos** que não pertence a um anel
            //3. Dado que (CH2=CH)-R entao o numero de elementos ligados em Ci e Cj tem que ser (Ci == 1 e Cj == 2) ou (Ci == 2 e Cj == 1)
            if(MX.mol[m].MAG[j][4] == 3 && MX.mol[m].MAG[j][2] == 999) 
            {
                if ((MX.mol[m].MAGnode[MX.mol[m].MAG[j][0]] == 2 ) && (MX.mol[m].MAGnode[MX.mol[m].MAG[j][1]] == 2 ))
                {
                    CH_A = MX.mol[m].MAG[j][0];
                    C = MX.mol[m].MAG[j][1];
                    // buscando as ligações deste átomo '=C='
                    //Liga_C = BuscaLigas_e(m, MX.mol[m].MAG[j][1]);
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C];
                    

                    for (unsigned int li = 0; li < Liga_C.size(); li++)
                    { 
                        if ( (Liga_C[li] > -1) && (MX.mol[m].MAG[Liga_C[li]][1] != MX.mol[m].MAG[j][1])) {

                            //verifica se e o Carbono (=C=) em ambas ligacoes e o mesmo, se for escolhe o seu companheiro da ligacao (CH).  
                            if (MX.mol[m].MAG[Liga_C[li]][0] == MX.mol[m].MAG[j][1] && (MX.mol[m].MAG[Liga_C[li]][2] == 999) && (MX.mol[m].MAG[Liga_C[li]][4] == 3) && (MX.mol[m].MAGnode[MX.mol[m].MAG[Liga_C[li]][1]] == 2))
                            { 
                                CH_B = MX.mol[m].MAG[Liga_C[li]][1];
                                if (resposta.size() == unsigned(1) && resposta[0][0] == -1 &&  resposta[0][1] == -1 && resposta[0][2] == -1)
                                {
                                    resposta[0][0] = CH_A;
                                    resposta[0][1] = C; 
                                    resposta[0][2] = CH_B;
                                } else {
                                    std::vector < int > _newMG = { CH_A, C, CH_B }; 
                                    resposta.push_back(std::move(_newMG));
                                }
                            } 
                        }
                    }  
                } 
            }
        }

        return resposta;
    }


    #pragma endregion Estruturas_1_ordem



    /////////////////////////////////////////////////////
    ////             FUNÇÕES AUXILIARES              ////
    /////////////////////////////////////////////////////
    #pragma region funcoes_auxiliares

    // Busca ligações que contÊm o elemento 'e' na molécula 'm'
    std::vector < int > BuscaLigas_e(int m, int e, Mistura &MX)
    {
        std::vector < int > resposta(4,-1);// {-1,-1,-1,-1};
        unsigned int n = 0;
        for (unsigned int L = 0 ; L < MX.mol[m].MAG.size() ; L++)
        {// varrendo todas as ligações 'L' da molécula 'm'
            if (MX.mol[m].MAG[L][0] == e || MX.mol[m].MAG[L][1] == e) 
            {
                resposta[n] = L;
                n++;
            }
            if (MX.mol[m].MAG[L][0] < 0) L = MX.mol[m].MAG.size();
        }
        return resposta;
    }

    // Testando se Pelo Menos Uma Ligação do Elemento 'e' da molécula 'm' é do tipo 't'
    bool pmul_e(int m, int e, int t, Mistura &MX)
    {
        // Variável de resposta do tipo. Se _tipo=1 elemento possui pelo menos uma ligação do tipo 't'. Se _tipo=0 não possui
        bool _tipo = 0;

        // Vetor que armazena o numero das ligações que o elemento 'e' participa
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // buscando as ligações deste elemento 'e'
        Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

        // varrendo ligações de 'e'
        for (unsigned int L = 0 ; L < Liga_e.size() ; L++)
        {
            // testando se esta ligação deste elemento 'e' é do tipo 't'
            if (MX.mol[m].MAG[Liga_e[L]][4] == t) _tipo = 1;
            else if (MX.mol[m].MAG[Liga_e[L]][4] >= 0) _tipo = 0;

            // se a ligação 'L' for no tipo 't' ou se não existe a próxima ligação termina o 'for'
            if ((Liga_e[L+1] < 0) || (_tipo == 1)) L = Liga_e.size();
        }
        return _tipo;
    }

    // Testando se Todas As Ligação do Elemento 'e' da molécula 'm' é do tipo 't'
    bool tal_e(int m, int e, int t, Mistura &MX)
    {
        // Variável de resposta do tipo. Se _tipo = 1 todas as ligações do elemento são do tipo 't'. Se _tipo = 0 ao menos uma não é do tipo
        bool _tipo = 0;

        // Vetor que armazena o numero das ligações que o elemento 'e' participa
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // buscando as ligações deste elemento 'e'
        Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

        // varrendo ligações de 'e'
        for (unsigned int L = 0 ; L < Liga_e.size() ; L++)
        {
            // testando se esta ligação deste elemento 'e' é do tipo 't'
            if (MX.mol[m].MAG[Liga_e[L]][4] == t) _tipo = 1;
            else if (MX.mol[m].MAG[Liga_e[L]][4] >= 0) _tipo = 0;

            // se a ligação 'L' não for no tipo 't' ou se não existe a próxima ligação termina o 'for'
            if ((Liga_e[L+1] < 0) || (_tipo == 0)) L = Liga_e.size();
        }
        return _tipo;
    }

    // Testando se Todas Os Elementos Ligados ao Elemento 'e' da molécula 'm' são Carbonos
    bool toel_e(int m, int e, Mistura &MX)
    {
        // Variável de resposta do tipo. Se _tipo = 1 todos os elementos ligados ao elemento 'e' são Carbonos. Se _tipo = 0 pelo menos 1 não é
        bool _tipo = 0;

        // Átomos que se ligam a E
        unsigned int A = 0;

        // Vetor que armazena o numero das ligações que o elemento 'e' participa
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // buscando as ligações deste elemento 'e'
        Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

        // varrendo ligações de 'e'
        for (unsigned int L = 0 ; L < Liga_e.size() ; L++)
        {
            // Pegando o primeiro elemento que se liga a 'e' 
            if (MX.mol[m].MAG[Liga_e[L]][0] == e) A = MX.mol[m].MAG[Liga_e[L]][1];
            else A = MX.mol[m].MAG[Liga_e[L]][0];

            // testando se é carbono
            if (MX.mol[m].MAGnode[A] < 1000) _tipo = 1;
            else _tipo = 0;

            // se A  não for carbono ou não houver uma próxima ligação interrompe po for
            if ((Liga_e[L+1] < 0) || (_tipo == 0)) L = Liga_e.size();
        }
        return _tipo;
    }

    // Conta o número de ligações do anel ligado
    int ContaLigaAnel(int m,int anel, Mistura &MX)
    {
        int contador = 0;

        //Varrendo o vetor de ligações do atomo A
        for (unsigned int i = 0; i < MX.mol[m].MAG.size();i++)
        {
            if (MX.mol[m].MAG[i][2] == anel) contador++;

            if (MX.mol[m].MAG[i+1][2] == -1) i = MX.mol[m].MAG.size();
        }
        return contador;
    }


    // Identificação do Anel predominante do Carbono
    int AnelPred(int m,int A, Mistura &MX)
    {
        int anel = 0;
        int Liga_anel = 0;
        int anel_max = 0;
        int Liga_anel_max = 0;
        std::vector < int > Liga_A (CTE::maxLigas_e, -1);


        // buscando as ligações deste átomo 'A'
        Liga_A = MX.mol[m].LigacoesPorAtomo[A];// BuscaLigas_e(m, A, MX);

        //Varrendo o vetor de ligações do atomo A
        for (unsigned int i = 0; i < Liga_A.size(); i++)
        {
            //Identificando o anel da ligação i
            anel = MX.mol[m].MAG[Liga_A[i]][2];
            Liga_anel = 0;

            //Conta as ligações do anel identificado acima 
            for (unsigned int j = 0; j < Liga_A.size(); j++) 
            {
                if (MX.mol[m].MAG[Liga_A[j]][2] == anel) Liga_anel++;
                
                //proteger o for quando A não tem todas as ligações, ie , apresenta -1
                if ((j < Liga_A.size()-1) && (Liga_A[j+1] == -1)) j = Liga_A.size();
            }
            
            // testa se tem mais ligações do que o máximo
            if (Liga_anel > Liga_anel_max) 
            {
                Liga_anel_max = Liga_anel;
                anel_max = anel;
            }

            if ((i < Liga_A.size()-1) && (Liga_A[i+1] == -1)) i = Liga_A.size();
        }
        return anel_max;
    }


    //determina o numero de aneis da molecula, em todos os nucleos.
    int nrRings_m(int m, Mistura &MX)
    {
        int resposta = -1;
        for (unsigned int L = 0 ; L < MX.mol[m].MAG.size(); L++)
        {
            if ( (MX.mol[m].MAG[L][2] < 999) && (MX.mol[m].MAG[L][3] < 998) && (MX.mol[m].MAG[L][2] > -1)) 
            {
                resposta = MX.mol[m].MAG[L][2]; 
                //if (MX.mol[m].MAG[L][2] > resposta) resposta += 1;
            } 
            else
            {
                continue;
            }

            if (MX.mol[m].MAG[L][0] < 0) L = MX.mol[m].MAG.size(); // Acabou a molécula

        }

        return resposta;
    }

    //determina o numero de aneis por nucleo.
    std::vector<int> nrRings_k_m(int m, Mistura &MX)
    {
        std::vector<int> resposta;
        resposta.resize(CTE::maxNucleos);
        for (auto &v:resposta) { v = 0; }
        std::vector<int> aneisContabilizados; 

        for (unsigned int L = 0 ; L < MX.mol[m].MAG.size(); L++) {
            if ( (MX.mol[m].MAG[L][2] < 999) && (MX.mol[m].MAG[L][3] < 998) && (MX.mol[m].MAG[L][2] > -1) 
            && std::find(aneisContabilizados.begin(), aneisContabilizados.end(), MX.mol[m].MAG[L][2]) == aneisContabilizados.end()) {
                aneisContabilizados.push_back(MX.mol[m].MAG[L][2]);
                resposta[MX.mol[m].MAG[L][3]] += 1; 
            } 
            else {
                continue;
            }

            if (MX.mol[m].MAG[L][0] < 0) L = MX.mol[m].MAG.size(); // Acabou a molécula

        }

        return resposta;
    }

    std::vector < std::vector < int > > BuscaLigas_Anel(int m, Mistura &MX)
    {
        int _nrRings = nrRings_m(m, MX);
        std::vector < std::vector < int > > resposta;
        resposta.resize(MX.mol[m].MAG.size());
        for (int a = 0 ; a < _nrRings + 1; a++)
        {
            resposta[a].resize(1);
            resposta[a][0] = -1;
            
            for (unsigned int L = 0 ; L < MX.mol[m].MAG.size() ; L++)
            {// varrendo todas as ligações 'L' da molécula 'm'
                if ( MX.mol[m].MAG[L][2] == a) 
                {
                    std::vector < int > Liga_C1(CTE::maxLigas_e, -1);

                    //grava todas as ligacoes do anel
                    if (resposta[a][0] == -1) {
                        resposta[a][0] = L;
                    } else {
                        resposta[a].push_back(L);
                    };

                    if ( resposta[a][0] != -1 )
                    {
                        //procura moleculas penduradas
                        if (MX.mol[m].MAGnode[MX.mol[m].MAG[L][0]] == 3 ) 
                        {
                            //CarbonosTerciarios.push_back(MX.mol[m].MAG[L][0]);
                            int _carbonoDoAnel = MX.mol[m].MAG[L][0]; //Ci  variavel local
                            //Liga_C1 = BuscaLigas_e(m, carbonoDoAnel, MX);
                            Liga_C1 = MX.mol[m].LigacoesPorAtomo[_carbonoDoAnel];
                            for (unsigned int iii = 0; iii < Liga_C1.size(); iii++)
                            {
                                //vai excluir se houver ligacao entre nucleos do tipo AA ? MX.mol[m].MAG[Liga_C1[iii]][3] == 998 ? 
                                if (Liga_C1[iii] != -1 && Liga_C1[iii] != signed(L) && MX.mol[m].MAG[Liga_C1[iii]][2] >= 999 
                                && ((MX.mol[m].MAG[Liga_C1[iii]][0] == _carbonoDoAnel) || (MX.mol[m].MAG[Liga_C1[iii]][1] == _carbonoDoAnel)))
                                { // escolher so posicao Ci ou Cj para registrar ? ou verificar se ja foi adicionado no vetor CarbonosTerciarios  ((MX.mol[m].MAG[Liga_C1[iii]][0] == carbonoDoAnel) || (MX.mol[m].MAG[Liga_C1[iii]][1] == carbonoDoAnel))) ?
                                    resposta[a].push_back(Liga_C1[iii]);
                                }
                            }
                        }

                    }

                } // caso da Piridina 
                else if ( (MX.mol[m].MAG[L][2] == 1001) && (pmul_e(m, MX.mol[m].MAG[L][0], 2, MX)) && (pmul_e(m, MX.mol[m].MAG[L][1], 2, MX))) 
                {
                    
                    if (resposta[a][0] == -1) {
                        resposta[a][0] = L;
                    } else {
                        resposta[a].push_back(L);
                    }
                }
                
                if (MX.mol[m].MAG[L][0] < 0) L = MX.mol[m].MAG.size(); //Acabou a molécula
            }
        }

        return resposta;

    }

    // numero de Ligação do Elemento 'e' da molécula 'm' é do tipo 't'
    int nrlt_e(int m, int e, int t, Mistura &MX)
    {
        // Variável de resposta do tipo. Se _tipo=1 elemento possui pelo menos uma ligação do tipo 't'. Se _tipo=0 não possui
        int _nrLigacoesTipo = 0;

        // Vetor que armazena o numero das ligações que o elemento 'e' participa
        std::vector < int > Liga_e (CTE::maxLigas_e, -1);

        // buscando as ligações deste elemento 'e'
        //Liga_e = BuscaLigas_e(m, e);
        Liga_e = MX.mol[m].LigacoesPorAtomo[e];//BuscaLigas_e(m, e, MX);

        // varrendo ligações de 'e'
        for (unsigned int L = 0 ; L < Liga_e.size() ; L++)
        {
            // testando se esta ligação deste elemento 'e' é do tipo 't'
            if (MX.mol[m].MAG[Liga_e[L]][4] == t) _nrLigacoesTipo += 1;

            // se a ligação 'L' for no tipo 't' ou se não existe a próxima ligação termina o 'for'
            if (Liga_e[L+1] < 0) L = Liga_e.size();
        }

        return _nrLigacoesTipo;
    }

    // Identificando se o CH2 pertence ao grupo (CH2-C=O - R)
    bool CH2pertenceKO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GCH2COCH2(m, MX);   

        // identificando se o CH2 está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][0] == A) _tipo = 1;
        }
        return _tipo;
    }

    // Identificando se o CH pertence ao grupo (CH-C=O - R)
    bool CHpertenceKO(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > GrupoCentona = GCHCO(m, MX);   

        // identificando se o CH está ligado a um grupo centona
        for (unsigned int i = 0; i < GrupoCentona.size(); i++)
        {
            if (GrupoCentona[i][0] == A) _tipo = 1;
        }
        return _tipo;
    }

    // Identificando se o aC pertence ao grupo aC-CH2-aC 
    bool aCpertence_aCCH2aC(int m, int A, Mistura &MX)
    {
        // variáveis auxiliares da função
        bool _tipo = 0;
        std::vector < std::vector < int > > Grupo = GaCCH2aC(m, MX);   

        // identificando se o CH2 está ligado a um grupo centona
        for (unsigned int i = 0; i < Grupo.size(); i++)
        {
            if (Grupo[i][0] == A) _tipo = 1;
        }
        return _tipo;
    }



    // Identificando o grupo aC-CH=O aldeido
    std::vector < std::vector < int > > Galdeido_aC_CO(int m, Mistura &MX)
    {
        int KO = 0;
        int cont = 0;
        int _tipo_aC = -1; 
        int _tipoNao_aC = -1;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);
            
        // Conta quantos grupos cetônicos(ou aldeido) tem na molécula 
        KO = MX.mol[m].SOLex[0][39];

        if (KO == 0) KO++;
        resposta.resize(KO);
        A.resize(KO);
        B.resize(KO);
        C.resize(KO);

        // criando a matriz de resposta do grupo aC-CH=O
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {   
            //encontrando o oxigenio do grupo centonico ou aldeido
            if(MX.mol[m].MAG[j][2] == 1003)
            {   
                //encontrando o carbono (C) ligado ao oxigenio do aldeído, e o carbonos ligado à "C", "A" ou "B".
                C[cont] = MX.mol[m].MAG[j][0];

                if (MX.mol[m].MAGnode[C[cont]] == 2)            
                {
                    // buscando as ligações deste átomo 'C'           
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
                        
                    if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                    else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                    if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                    else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];

                    _tipoNao_aC = -1;
                    // testando se o carbono "A" é terciário e aromático e a outra ligação tem que ser com o heteroátomo
                    if(MX.mol[m].MAGnode[A[cont]] == 3 && pmul_e(m, A[cont], 2, MX) == 1 && MX.mol[m].MAGnode[B[cont]] == 1003)                
                    {
                        _tipo_aC = A[cont]; 
                        _tipoNao_aC = B[cont];

                    }    
                    // testando se o carbono "B" é terciário e aromático e a outra ligação tem que ser com o heteroátomo.
                    else if (MX.mol[m].MAGnode[B[cont]] == 3 && pmul_e(m, B[cont], 2, MX) == 1 && MX.mol[m].MAGnode[A[cont]] == 1003)
                    {
                        _tipo_aC = B[cont];
                        _tipoNao_aC = A[cont];
                    }

                    //O carbono R,  tem que existir para que exista o grupo aC-CO-R.
                    if (_tipoNao_aC >= 0 )                
                    {
                        resposta[cont][0] = _tipo_aC;
                        resposta[cont][1] = C[cont];
                        resposta[cont][2] = _tipoNao_aC;
                        cont++;
                    } 
                }       
            } 
        }
        return resposta;
    }

    // Identificação a estrutura da cetona -->  CH-(C=O)-R
    std::vector < std::vector < int > > GCHCO(int m, Mistura &MX)
    {
        int KO = 0;
        int cont = 0;
        int _tipoCH = -1; 
        int _tipoNaoCH = -1;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);

        //Conta quantos grupos cetônicos tem na molécula 
        KO = MX.mol[m].SOLex[0][39];

        if (KO == 0) KO++;
        resposta.resize(KO);
        A.resize(KO);
        B.resize(KO);
        C.resize(KO);

        //criando a matriz de resposta do grupo CH-(C=O)-R
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {   
            //encontrando o oxigenio do grupo centonico
            if(MX.mol[m].MAG[j][2] == 1003)
            {   
                //encontrando o carbono (C) ligado ao oxigenio da cetona, e os carbonos ligados à "C", "A" e "B".
                C[cont] = MX.mol[m].MAG[j][0];

                if (MX.mol[m].MAGnode[C[cont]] == 3)            
                {
                    // buscando as ligações deste átomo 'C'           
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
                        
                    if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                    else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                    if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                    else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];

                    _tipoNaoCH = -1;
                    // testando se o carbono "A" é terciário e não aromático.
                    if(MX.mol[m].MAGnode[A[cont]] == 3 && pmul_e(m, A[cont], 2, MX) == 0 && MX.mol[m].MAGnode[B[cont]] < 1000)
                    {
                        _tipoCH = A[cont]; 
                        _tipoNaoCH = B[cont];

                    }    
                    // testando se o carbono "B" é terciário e não aromático.
                    else if (MX.mol[m].MAGnode[B[cont]] == 3 && pmul_e(m, B[cont], 2, MX) == 0 && MX.mol[m].MAGnode[A[cont]] < 1000)
                    {
                        _tipoCH = B[cont];
                        _tipoNaoCH = A[cont];
                    }

                    //O carbono que não é primário tem que existir e não ser aromático para construir a matriz dos grupos CH-CO-R
                    if (_tipoNaoCH >= 0 && pmul_e(m, _tipoNaoCH, 2, MX) == 0 && MX.mol[m].MAGnode[_tipoNaoCH] > 1)                
                    {
                        resposta[cont][0] = _tipoCH;
                        resposta[cont][1] = C[cont];
                        resposta[cont][2] = _tipoNaoCH;
                        cont++;
                    } 
                }       
            } 
        }
        return resposta;
    }

    // Identificação a estrutura desse tipo cetona--> CH2-(C=O)-CH2
    std::vector < std::vector < int > > GCH2COCH2(int m, Mistura &MX)
    {
        int KO = 0;
        int cont = 0;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);
        std::vector < std::vector < int > > resposta;

        //Conta quantos grupos cetônicos tem na molécula 
        KO = MX.mol[m].SOLex[0][39];

        if (KO == 0) KO++;
        resposta.resize(KO);
        A.resize(KO);
        B.resize(KO);
        C.resize(KO);

        //criando a matriz de resposta do grupo CH2-(C=O)-CH2
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {
            if(MX.mol[m].MAG[j][2] == 1003)
            {
                C[cont] = MX.mol[m].MAG[j][0];
                
                if (MX.mol[m].MAGnode[C[cont]] == 3)    
                {    
                    // buscando as ligações deste átomo 'C'           
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
                        
                    if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                    else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                    if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                    else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];
                
                    if(MX.mol[m].MAGnode[A[cont]] == 2 && MX.mol[m].MAGnode[B[cont]] == 2)
                    {
                        resposta[cont][0] = A[cont];
                        resposta[cont][1] = C[cont];
                        resposta[cont][2] = B[cont];
                        cont++;
                    } 
                }
            
            } 
        }
        return resposta;
    }

    // Identificação a estrutura da cetona -->  CH3-(C=O)-R
    std::vector < std::vector < int > > GCH3CO(int m, Mistura &MX)
    {
        int KO = 0;
        int cont = 0;
        int _tipoCH3 = -1; 
        int _tipoNaoCH3 = -1;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);
            
        //Conta quantos grupos cetônicos tem na molécula 
        KO = MX.mol[m].SOLex[0][39];

        if (KO == 0) KO++;
        resposta.resize(KO);
        A.resize(KO);
        B.resize(KO);
        C.resize(KO);

        //criando a matriz de resposta do grupo CH3-(C=O)-R
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {   
            //encontrando o oxigenio do grupo centonico
            if(MX.mol[m].MAG[j][2] == 1003)
            {   
                //encontrando o carbono (C) ligado ao oxigenio da cetona, e os carbonos ligados à "C", "A" e "B".
                C[cont] = MX.mol[m].MAG[j][0];

                if (MX.mol[m].MAGnode[C[cont]] == 3)            
                {
                    // buscando as ligações deste átomo 'C'           
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
                        
                    if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                    else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                    if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                    else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];

                    _tipoNaoCH3 = -1;
                    // testando se o carbono "A" é primário e o carbono "C" da cetona tem que ser terciário.
                    if(MX.mol[m].MAGnode[A[cont]] == 1 && MX.mol[m].MAGnode[B[cont]] < 1000)
                    {
                        _tipoCH3 = A[cont]; 
                        _tipoNaoCH3 = B[cont];
                    }    
                    // testando se o carbono "B" é primário e o carbono "C" da cetona tem que ser terciário.
                    else if (MX.mol[m].MAGnode[B[cont]] == 1 && MX.mol[m].MAGnode[A[cont]] < 1000) 
                    {
                        _tipoCH3 = B[cont];
                        _tipoNaoCH3 = A[cont];
                    }

                    //O carbono que não é primário tem que existir e não ser aromático para construir a matriz dos grupos CH3-CO-R
                    if (_tipoNaoCH3 >= 0 && pmul_e(m, _tipoNaoCH3, 2, MX) == 0)                
                    {
                        resposta[cont][0] = _tipoCH3;
                        resposta[cont][1] = C[cont];
                        resposta[cont][2] = _tipoNaoCH3;
                        cont++;
                    } 
                }
            } 
        }
        return resposta;
    }

    // Identificação a estrutura da cetona -->  aC-(C=O)-R
    std::vector < std::vector < int > > GaC_CO(int m, Mistura &MX)
    {
        int KO = 0;
        int cont = 0;
        int _tipo_aC = -1; 
        int _tipoNao_aC = -1;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);
            
        //Conta quantos grupos cetônicos tem na molécula 
        KO = MX.mol[m].SOLex[0][39];

        if (KO == 0) KO++;
        resposta.resize(KO);
        A.resize(KO);
        B.resize(KO);
        C.resize(KO);

        //criando a matriz de resposta do grupo aC-(C=O)-R
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAG.size(); j++) 
        {   
            //encontrando o oxigenio do grupo centonico
            if(MX.mol[m].MAG[j][2] == 1003)
            {   
                //encontrando o carbono (C) ligado ao oxigenio da cetona, e os carbonos ligados à "C", "A" e "B".
                C[cont] = MX.mol[m].MAG[j][0];

                if (MX.mol[m].MAGnode[C[cont]] == 3)            
                {
                    // buscando as ligações deste átomo 'C'           
                    Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
                        
                    if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                    else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                    if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                    else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];

                    _tipoNao_aC = -1;
                    // testando se o carbono "A" é terciário e aromático
                    if(MX.mol[m].MAGnode[A[cont]] == 3 && pmul_e(m, A[cont], 2, MX) == 1 && MX.mol[m].MAGnode[B[cont]] < 1000)                
                    {
                        _tipo_aC = A[cont]; 
                        _tipoNao_aC = B[cont];

                    }    
                    // testando se o carbono "B" é terciário .
                    else if (MX.mol[m].MAGnode[B[cont]] == 3 && pmul_e(m, B[cont], 2, MX) == 1 && MX.mol[m].MAGnode[A[cont]] < 1000)
                    {
                        _tipo_aC = B[cont];
                        _tipoNao_aC = A[cont];
                    }

                    //O carbono R,  tem que existir para que exista o grupo aC-CO-R.
                    if (_tipoNao_aC >= 0 )                
                    {
                        resposta[cont][0] = _tipo_aC;
                        resposta[cont][1] = C[cont];
                        resposta[cont][2] = _tipoNao_aC;
                        cont++;
                    } 
                }       
            } 
        }
        return resposta;
    }

    // Identifica o grupo aC-CH2-aC
    std::vector < std::vector < int > > GaCCH2aC(int m, Mistura &MX)
    {
        int cont = 0;
        int _tipo_aC1 = -1; 
        int _tipo_aC2 = -1;
        std::vector < int > A ;
        std::vector < int > B ;   
        std::vector < int > C ;
        std::vector < std::vector < int > > resposta;
        std::vector < int > Liga_C (CTE::maxLigas_e, -1);
            
        resposta.resize(10);
        A.resize(10);
        B.resize(10);
        C.resize(10);

        //criando a matriz de resposta do grupo aC-(C=O)-R
        for (unsigned int i = 0; i < resposta.size(); i++) 
        {
            resposta[i].resize(3);
            resposta[i][0] = -1;
            resposta[i][1] = -1;
            resposta[i][2] = -1;
        }
        // varrendo as linhas da matriz ligas
        for (unsigned int j = 0; j < MX.mol[m].MAGnode.size(); j++) 
        {   
            //encontrando o carbono CH2
            if(MX.mol[m].MAGnode[j] == 2 && tal_e(m, j, 0, MX))
            {   
                //encontrando o carbono (C) ligado ao oxigenio da cetona, e os carbonos ligados à "C", "A" e "B".
                // buscando as ligações deste átomo 'C' 
                C[cont] = j;
            
                Liga_C = MX.mol[m].LigacoesPorAtomo[C[cont]];//BuscaLigas_e(m, C[cont], MX); 
            
                                    
                if (MX.mol[m].MAG[Liga_C[0]][0] ==  C[cont]) A[cont] = MX.mol[m].MAG[Liga_C[0]][1];
                else  A[cont] = MX.mol[m].MAG[Liga_C[0]][0];

                if (MX.mol[m].MAG[Liga_C[1]][0] ==  C[cont]) B[cont] = MX.mol[m].MAG[Liga_C[1]][1];
                else  B[cont] = MX.mol[m].MAG[Liga_C[1]][0];

                
                // testando se o carbono "A" e B" têm que ser aromáticos
                if(pmul_e(m, A[cont], 2, MX) == 1 && pmul_e(m, B[cont], 2, MX) == 1)                
                {
                    _tipo_aC1 = A[cont]; 
                    _tipo_aC2 = B[cont];

                    resposta[cont][0] = _tipo_aC1;
                    resposta[cont][1] = C[cont];
                    resposta[cont][2] = _tipo_aC2;
                    cont++;
                } 
            } 
        }
        return resposta;
    }



    #pragma endregion funcoes_auxiliares

}

