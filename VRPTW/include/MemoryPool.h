#include <list>
#include "safe_vector.h"
#include <sys/mman.h>
#include <unistd.h>
#include <cassert>
#include <cstring>
#include <memory>

//#define assertm(exp, msg) assert(((void)msg, exp))
#define assertTrue(msg) assert(((void)msg, false))


namespace MemoryPool_NS
{
    template <typename T>
    using IteratorListT = std::list<T*>::iterator;

    template<typename T>
    using referenceT = std::reference_wrapper<T>;

#define DefRefW_T(T) using RefW_##T = std::reference_wrapper<T>;

    template <typename T>
    class ProxT_notUsed
    {
    public:
        IteratorListT<T> iterator;
        long proxPos = -1;
    };



    template <typename T>
    class Pool
    {
    private:
        std::list<T*> listT;
        ProxT_notUsed<T> proxT;                            // Aponta para um vetor e a pos da prox sol nao utilizada

        T **p_tDelT           = nullptr;                   // Solucoes que nao sao mais utilizadas.
        long numVetSolDel     = 0;
        long nextPosVetSolDel = 0;
        long pageSize         = 0;
        long bucketSize       = 0;
        long numElemBucket    = 0;
        bool poolStart        = false;

        void rmPool()
        {

            for(auto t:listT)
                free(t);

            free(p_tDelT);

            listT.erase(listT.begin(), listT.end());
            p_tDelT = nullptr;
        }

    public:

        Pool() = default;
        Pool(Pool<T> &p)                  = delete;
        Pool(Pool<T> &&p)                 = delete;
        void operator=(Pool<T> &p)        = delete;
        void operator=(Pool<T> &&p)       = delete;
        Pool(const Pool<T> &p)            = delete;
        Pool(const Pool<T> &&p)           = delete;
        void operator=(const Pool<T> &p)  = delete;
        void operator=(const Pool<T> &&p) = delete;
        ~Pool(){rmPool();}
        long getNumElemBucket()const{return numElemBucket;}
        long getNumDel()const{return nextPosVetSolDel;}
        long getNumBuckets()const{return listT.size();}

        Pool(int bucketSizePage, int vetDelSizePage)
        {
            startPool(bucketSizePage, vetDelSizePage);
        }

        void startPool(const int bucketSizePage=1, const int vetDelSizePage=4)
        {
            if(poolStart)
                return;

            if(bucketSizePage < 1 || vetDelSizePage < 1)
                throw std::bad_alloc();


            std::cout<<"Pool()\n";
            pageSize = sysconf(_SC_PAGE_SIZE);
            bucketSize = bucketSizePage*pageSize;

            std::cout<<"pageSize: "<<pageSize<<"\n";
            std::cout<<"bucketSize: "<<bucketSize<<"\n";

            T *p_t = nullptr;
            assertm(posix_memalign((void**)&p_t, pageSize, bucketSize), "posix_memalign");
            assertm(posix_memalign((void**)&p_tDelT, pageSize, vetDelSizePage*bucketSize), "posix_memalign");



//#if VAR_POOL_SET_MEM_PAGES == 1
            std::cout<<"memset pages!\n";
            memset(p_t, 0, bucketSize);
            memset(p_tDelT, 0, vetDelSizePage*bucketSize);
//#endif

            listT.push_back(p_t);



            std::printf("p_t: %p\n", p_t);
            //std::printf("p_tLast: %p\n", p_tLast);


            numElemBucket  = bucketSize/sizeof(T);
            numVetSolDel   = (vetDelSizePage*bucketSize)/sizeof(T*);

            std::cout<<"numElemBucket: \t"<<numElemBucket<<"\n";
            std::cout<<"numVetSolDel: \t"<<numVetSolDel<<"\n";

            proxT.iterator = listT.begin();
            proxT.proxPos  = 0;
            poolStart = true;
        }


        [[nodiscard]]
        T* getT()
        {
            assertm(!poolStart, "Pool didn't start!");

            T* p_tTemp = nullptr;

            const long lastPosDel = nextPosVetSolDel-1;
            if(lastPosDel >= 0)
            {
                //std::cout<<"From p_tDelT\n";
                p_tTemp = p_tDelT[lastPosDel];
                p_tDelT[lastPosDel] = nullptr;
                nextPosVetSolDel = lastPosDel;
                return p_tTemp;
                //return p_tTemp;
            }

            if(proxT.iterator == listT.end())
            {

                //std::cout<<"Net bucket\n";
                T *p_tNew = nullptr;
                assertm(posix_memalign((void**)&p_tNew, pageSize, bucketSize), "Bay more memory!");

//#ifdef POOL_SET_MEM_PAGES
                memset(p_tNew, 0, bucketSize);
//#endif

                listT.push_back(p_tNew);
                proxT.iterator = (--listT.end());
                proxT.proxPos = 0;
            }

            if(proxT.proxPos < numElemBucket)
            {
                p_tTemp = ((*proxT.iterator) + proxT.proxPos);
                proxT.proxPos += 1;

                if(proxT.proxPos == numElemBucket)
                {
                    proxT.proxPos = 0;
                    ++proxT.iterator;
                }

                return p_tTemp;

            }

            assertTrue("Erro, nao deveria chegar aqui!");
        }

        void resetPool(bool delPool)
        {
            assertm(!poolStart, "Pool didn't start!");

            if(delPool)
            {
                rmPool();
                poolStart = false;
            }

            proxT.iterator = listT.begin();
            proxT.proxPos  = 0;
            nextPosVetSolDel = 0;
        }

        void delT(T* p_t)
        {
            assertm(!poolStart, "Pool didn't start!");

            //std::cout<<"Del: "<<*p_t<<"\n";
            if(nextPosVetSolDel < numVetSolDel)
            {
                p_tDelT[nextPosVetSolDel] = p_t;
                nextPosVetSolDel += 1;
            }

            //p_tDelT[nextPosVetSolDel] = p_t;
            //nextPosVetSolDel += (nextPosVetSolDel < numVetSolDel);

        }
    };

}
