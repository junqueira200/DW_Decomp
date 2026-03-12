module;


#include <iostream>

export module teste;

//export namespace teste
//{
    export template <typename T>
    void func()
    {
        T val;
        std::cout<<"vall: "<<val<<"\n";
    }
//}