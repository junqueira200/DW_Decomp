//
// Created by igor on 19/09/24.
//

#ifndef DW_DECOMP_AUX_H
#define DW_DECOMP_AUX_H


#define __PRETTYFILE__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define PRINT_DEBUG(inicio, texto) std::cout<<inicio<<"DEBUG: "<<texto<<"  FILE: "<<__PRETTYFILE__<<"  FUNC: "<<__PRETTY_FUNCTION__<<"  LINHA: "<<__LINE__<<"\n";

#define assertm(exp, msg) if(exp){std::cout<<msg<<"\n\nFILE: "<<__PRETTYFILE__<<"\nFUNC: "<<__PRETTY_FUNCTION__<<"\n\n"; throw "ERROR";}

#endif //DW_DECOMP_AUX_H
