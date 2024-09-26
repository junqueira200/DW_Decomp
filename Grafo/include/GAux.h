//
// Created by igor on 25/09/24.
//

#ifndef MNFP_GAUX_H
#define MNFP_GAUX_H


#ifndef DEBUG_AUX
#define DEBUG_AUX
#define __PRETTYFILE__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define PRINT_DEBUG(inicio, texto) std::cout<<inicio<<"DEBUG: "<<texto<<"  FILE: "<<__PRETTYFILE__<<"  FUNC: "<<__PRETTY_FUNCTION__<<"  LINHA: "<<__LINE__<<"\n";


#endif // DEBUG_AUX
#endif //MNFP_GAUX_H
