/* ****************************************
 * ****************************************
 *  Data:    21/02/26
 *  Arquivo: MILP.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L_SDVRP
 * ****************************************
 * ****************************************/

#ifndef MILP_H
#define MILP_H

#include "gurobi_c++.h"
#include "AuxT.h"
#include "Solucao.h"


namespace MILP_NS
{

    class VectorGRBVar
    {
    private:
        bool    started         = false;
        GRBVar *vetVar          = nullptr;
        double *vetDoubleAttr_X = nullptr;
        int     num             = 0;
        char    varType         = GRB_CONTINUOUS;

    public:

        inline __attribute__((always_inline)) int getNum(){return num;}
        VectorGRBVar(GRBModel &model, int num_, const std::string &&name, char type);
        void start(GRBModel &model, int num_, const std::string &&name, char type);

        VectorGRBVar()=default;
        void operator = (const VectorGRBVar &)=delete;
        ~VectorGRBVar(){delete []vetVar; delete []vetDoubleAttr_X;}
        inline __attribute__((always_inline)) GRBVar& operator ()(const int indexI)
        {

        #if VAR_VECTOR_SANITY_CHECK
            if(indexI >= num || indexI < 0)
            {
                std::cout<<"Error indice i: "<<indexI<<" is wrong to vector of size: "<<num<<"\n";
                throw std::out_of_range("");
            }
        #endif

            return vetVar[indexI];
        }

        inline __attribute__((always_inline)) double getX_value(int indexI)
        {

            if(vetDoubleAttr_X == nullptr)
            {
                std::cout<<"ERRO, vetDoubleAttr is equal to nullptr\n";
                PRINT_DEBUGG("", "");
                exit(-1);
            }

        #if VAR_VECTOR_SANITY_CHECK
            if(indexI >= num || indexI < 0)
            {
                std::cout<<"Error indice i: "<<indexI<<" is wrong to vector of size: "<<num<<"\n";
                throw std::out_of_range("");
            }
        #endif

            return vetDoubleAttr_X[indexI];
        }

        void setUB(double ub);
        void setLB(double lb);
        void setUB_LB(double ub, double lb);
        void setUB_LB(double ub, double lb, int i);
        void printVars();

        void setVetDoubleAttr_X(GRBModel &model, bool Xn);
        void setAttr_Start0();
    };

    class MatrixGRBVar
    {
    private:
        bool inicializado = false;
        GRBVar *vetVar = nullptr;
        double *vetDoubleAttr_X = nullptr;
        int numLin, numCol;
        char typeVar;

    public:

        MatrixGRBVar(GRBModel &model, int numLin, int numCol, const std::string &&name, char type, bool zeroI_EqualJ);
        MatrixGRBVar()=default;
        MatrixGRBVar(const MatrixGRBVar&)=delete;
        void operator = (MatrixGRBVar &)=delete;
        void start(GRBModel &model, int numLin, int numCol, const std::string &&nome, char type, bool zeroI_EqualJ);

        ~MatrixGRBVar(){delete []vetVar;};
        inline __attribute__((always_inline)) GRBVar& operator ()(const int indexI, const int indexJ)
        {

#if VAR_VECTOR_SANITY_CHECK
            if(indexI >= numLin || indexI < 0)
            {
                std::cout<<"Erro indice i: "<<indexI<<" esta errado para matrix de tam "<<numLin<<" x "<<numCol<<"\n";
                throw std::out_of_range("");
            }

            if(indexJ >= numCol)
            {
                std::cout<<"Erro indice j: "<<indexJ<<" esta errado para matrix de tam "<<numLin<<" x "<<numCol<<"\n";
                throw std::out_of_range("");
            }
#endif

            return vetVar[indexI*numCol+indexJ];
        }

        inline __attribute__((always_inline)) double getX_value(const int indexI, const int indexJ)
        {
            if(vetDoubleAttr_X == nullptr)
            {
                std::cout<<"ERRO, vetDoubleAttr eh igual a nullptr\n";

            }


#if VAR_VECTOR_SANITY_CHECK
            if(indexI >= numLin || indexI < 0)
            {
                std::cout<<"Erro indice i: "<<indexI<<" esta errado para matrix de tam "<<numLin<<" x "<<numCol<<"\n";
                throw std::out_of_range("");
            }

            if(indexJ >= numCol)
            {
                std::cout<<"Erro indice j: "<<indexJ<<" esta errado para matrix de tam "<<numLin<<" x "<<numCol<<"\n";
                throw std::out_of_range("");
            }
#endif

            return vetDoubleAttr_X[indexI * numCol + indexJ];
        }

        void setUB(double ub);
        void setLB(double lb);
        void setUB_LB(double ub, double lb);
        void printVars();

        void setVetDoubleAttr_X(GRBModel &model, bool X_n);
        void setAttr_Start0();

    };

    class Variables
    {
    public:

        VectorGRBVar vetPosX;
        VectorGRBVar vetPosY;
        VectorGRBVar vetPosZ;

        VectorGRBVar vetDX;
        VectorGRBVar vetDY;
        VectorGRBVar vetDZ;

        MatrixGRBVar matX_pos;   // indicates if item i is placed to the left of item j
        MatrixGRBVar matY_pos;   // indicates if item i is placed in behind of item j
        MatrixGRBVar matZ_pos;   // indicates if item i is placed below of item j

        MatrixGRBVar matX_neg;   // indicates if item i is placed to the right of item j
        MatrixGRBVar matY_neg;   // indicates if item i is placed front of item j
        MatrixGRBVar matZ_neg;   // indicates if item i is placed above of item j


        MatrixGRBVar matRot;    // Indicates for the item i if the j rotation is used


        Variables(GRBModel& model, VectorI& vetItems, int numItems);

    };

    void addBasicConstraints(GRBModel& model, Variables& variables, SolucaoNS::Bin& bin);


}

#endif // MILP_H
