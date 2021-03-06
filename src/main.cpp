#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstring>
#include <list>
#include <set>
#include <algorithm>

using namespace std;


/************************ Global Variables ***************************/

// size
int n;
// fow matrix
vector< vector <int> > f;
// distances matrix
vector< vector <int> > d;


/****************************** Utils Functions *********************************/

float Rand(void){
    return ((double) rand() / (RAND_MAX));
}

int Randint(int low, int high){
    return rand() % (low+high) + low;
}

inline void swap(int &a, int &b){
    int aux = a;
    a = b;
    b = aux;
}


/****************************** Sort struct *********************************/

// Sort an array of pair in decreasing orden by second element
struct bySecond1{
    bool operator () (const pair<int,int>& a, const pair<int,int>& b){
       return a.second >=b.second;
    }
};

// Sort an array of pair in increasing orden by second element
struct bySecond2{
    bool operator () (const pair<int,int>& a, const pair<int,int>& b){
       return a.second <=b.second;
    }
};


/****************************** Cost Function *********************************/

// Calculate a cost
int costFunction(const vector<int> &s){
    int sum = 0;

    for(int j = 0; j < n; j++){
        for(int k = 0; k < n; k++){
            sum += f[j][k]*d[s[j]][s[k]];
        }
    }

    return sum;
}

// Calculate the diference between two solutions which differ in two positions i, j
int costDif(const vector<int> &s, int i, int j){
    int sum = 0;

    for(int k = 0; k < n; k++){
        if(k != i && k != j){
            sum+= f[k][i] * ( d[s[k]][s[j]] - d[s[k]][s[i]] ) +
                  f[k][j] * ( d[s[k]][s[i]] - d[s[k]][s[j]] ) +
                  f[i][k] * ( d[s[j]][s[k]] - d[s[i]][s[k]] ) +
                  f[j][k] * ( d[s[i]][s[k]] - d[s[j]][s[k]] );
        }
    }

    return sum;
}


/****************************** GREDDY *********************************/

// GREDDY SOLUTION
vector<int> greddy_QAP(){

    set<pair<int,int>, bySecond1> pf;
    set<pair<int,int>, bySecond2> pd;

    // Calculate potential flow (pf) and potential distance (pd)
    for(int i = 0; i < n; i++){
        int sumf = 0;
        int sumd = 0;
        for(int j = 0; j < n; j++){
            sumf += f[i][j];
            sumd += d[i][j];
        }
        pf.insert(pair<int,int>(i, sumf));
        pd.insert(pair<int,int>(i, sumd));
    }

    vector<int> s(n);

    for(int i = 0; i < n; i++){
        s[(pf.begin())->first] = (pd.begin())->first;
        pf.erase(pf.begin());
        pd.erase(pd.begin());
    }

    return s;
}


/****************************** LOCAL SEARCH *********************************/

// BL-QAP
vector<int> bl_QAP(){
    vector<int> s;
    vector<int> dlb(n, false);
    bool improve_flag = true;
    int iters = 0;

    for(int i = 0; i < n; i++)
        s.push_back(i);

    int cost = costFunction(s);

    while(improve_flag && iters < 50000){
        for(int i = 0; i < n; i++){
            if(dlb[i] == 0){
                improve_flag = false;
                for(int j = 0; j < n; j++){
                    int cd = costDif(s, i, j);
                    if(cd < 0){
                        swap(s[i], s[j]);
                        cost += cd;
                        dlb[i] = 0;
                        dlb[j] = 0;
                        improve_flag = true;
                    }
                    iters++;
                }
                if(!improve_flag)
                    dlb[i] = 1;
            }
        }
    }

    return s;
}


/********************************** GENETIC UTILS ************************************/

struct chromosome{
    // Permutation vector
    vector<int> s;
    // Fitness value
    double value;

    chromosome(vector<int> _s, double v){
        s = _s;
        value = v;
    }
    chromosome(){}
};

bool comp_chromosome(const chromosome &c1, const chromosome &c2){
    return c1.value < c2.value;
}


chromosome crossChromosome(const chromosome &c1, const chromosome &c2){

    int s_size = c1.s.size();

    vector<int> c(s_size, -1);
    vector<int> restos;

    for(int i = 0; i < s_size; i++){
        if(c1.s[i] == c2.s[i]){
            c[i] = c1.s[i];
        }
        else{
            restos.push_back(c1.s[i]);
        }
    }

    int restos_i = 0;

    for(int i = 0; i < s_size; i++){
        if(c[i] == -1){
            c[i] = restos[restos_i];
            restos_i++;
        }
    }

    return chromosome(c, costFunction(c));
}


void competition(vector<chromosome> &population, chromosome &child1, chromosome & child2){

    int tam_p = population.size();

    sort(population.begin(), population.end(), comp_chromosome);

    if(child1.value >= child2.value){
        if(population[tam_p-1].value > child1.value){
            population[tam_p-1] = child1;
            if(population[tam_p-2].value > child2.value){
                population[tam_p-2] = child2;
            }
        }
    }
    else if(population[tam_p-1].value > child2.value){
        population[tam_p-1] = child2;
        if(population[tam_p-2].value > child1.value){
            population[tam_p-2] = child1;
        }
    }

}

chromosome binaryTournament(const vector<chromosome> &population){
    chromosome father;
    int k1, k2;
    int tam_p = population.size();

    k1 = Randint(0, tam_p-1);
    do{k2 = Randint(0, tam_p-1);}while(k1 == k2);

    if(population[k1].value < population[k2].value)
        father = population[k1];
    else
        father = population[k2];


    return father;
}

void initPopulation(vector<chromosome> &population, int tam_p, int n){
    vector<int> s;

    for(int i = 0; i < n; i++){
        s.push_back(i);
    }

    for(int i = 0; i < tam_p; i++){
        random_shuffle(s.begin(), s.end());
        population.push_back(chromosome(s, costFunction(s)));
    }
}


/********************************** STATIONARY ************************************/

vector<int> Stationary(int max_generations = 1000){
    int tam_p = 50;
    vector<chromosome> population;
    initPopulation(population, tam_p, n);

    for(int generation = 0; generation < max_generations; generation++){
        chromosome father1, father2, child1, child2;

        father1 = binaryTournament(population);
        father2 = binaryTournament(population);
        child1 = crossChromosome(father1, father2);

        father1 = binaryTournament(population);
        father2 = binaryTournament(population);
        child2 = crossChromosome(father1, father2);

        for(int i = 0; i < n; i++){
            if(Rand() < 0.001){
                int rand_index = Randint(0, n);
                swap(child1.s[i], child1.s[rand_index]);
                child1.value = costFunction(child1.s);
            }
            if(Rand() < 0.001){
                int rand_index = Randint(0, n);
                swap(child2.s[i], child2.s[rand_index]);
                child2.value = costFunction(child2.s);
            }
        }
        competition(population, child1, child2);
    }

    int best = 0;

    for(int i = 1; i < tam_p; i++)
        if(population[i].value < population[best].value)
            best = i;

    return population[best].s;
}

/********************************** BL CHROMOSOME ************************************/

// BL-QAP
chromosome bl_QAP(chromosome c){
    vector<int> s = c.s;
    vector<int> dlb(n, false);
    bool improve_flag = true;
    int iters = 0;

    int cost = c.value;

    while(improve_flag && iters < 100){
        for(int i = 0; i < n; i++){
            if(dlb[i] == 0){
                improve_flag = false;
                for(int j = 0; j < n; j++){
                    int cd = costDif(s, i, j);
                    if(cd < 0){
                        swap(s[i], s[j]);
                        cost += cd;
                        dlb[i] = 0;
                        dlb[j] = 0;
                        improve_flag = true;
                    }
                    iters++;
                }
                if(!improve_flag)
                    dlb[i] = 1;
            }
        }
    }

    return chromosome(s, costFunction(s));
}

/********************************** GENERATIONAL ************************************/

vector<int> Generational(string algoritmo, int max_generations = 100){
    // Population size
    int tam_p = 50;
    // Número máximo de generaciones
    vector<chromosome> population;

    int best = 0, worst = 0;

    // init population
    initPopulation(population, tam_p, n);

    // Buscamos al mejor y peor de la población
    for(int i = 1; i < tam_p; i++){
        if(population[i].value < population[best].value){
            best = i;
        }
        if(population[i].value > population[worst].value){
            worst = i;
        }
    }

    // Iterate generations
    for(int generations = 1; generations < max_generations; generations++){
        // Save the best individual in the population
        chromosome best_chromosome_old = population[best];

        // Aplicar variantes con búsqueda local
        if(algoritmo != "ESTANDAR" && generations%10 == 0){
            // Algoritmo baldwiniano: solo obtenemos el fitness ya que el
            // cruce se hace sin las características aprendidas
            if(algoritmo == "BALDWINIANO"){
                for(int j = 0; j < tam_p; j++){
                    population[j].value = bl_QAP(population[j]).value;
                }
            }
            // Algoritmo lamarckiano: esta vez los individuos si que transmiten
            // las características aprendidas
            else if(algoritmo == "LAMARCKIANO"){
                for(int j = 0; j < tam_p; j++){
                    population[j] = bl_QAP(population[j]);
                }
            }
        }

        for(int j = 0; j < tam_p; j++){
            // 0.7 probability of cross
            if(Rand() < 0.7){
                chromosome father1, father2, child;

                // Get parents by bnary tournament
                father1 = binaryTournament(population);
                father2 = binaryTournament(population);

                // Cross parents to get the child
                child = crossChromosome(father1, father2);

                for(int k = 0; k < n; k++){
                    // 0.001 probability of gen mutation
                    if(Rand() <= 0.001){
                        int rand_index = Randint(0, n);
                        swap(child.s[k], child.s[rand_index]);
                        child.value = costFunction(child.s);
                    }
                }
                population[j] = child;
            }

            // Get the index of best individual in the population
            if(population[j].value < population[best].value){
                best = j;
            }
            // Get the index of worst individual in the population
            if(population[j].value > population[worst].value){
                worst = j;
            }
        }
        // Replace the worst of the population
        population[worst] = best_chromosome_old;

        if(best_chromosome_old.value < population[best].value)
            best = worst;
    }

    return bl_QAP(population[best]).s;
}

/********************************************************************************/


void input(char* nombre_fichero){
    ifstream fe(nombre_fichero);
    fe >> n;

    // Introduce flow
    for(int i = 0; i < n; i++){
        vector<int> v(n);
        for(int j = 0; j < n; j++)
            fe >> v[j];
        f.push_back(v);
    }

    // Introduce distances
    for(int i = 0; i < n; i++){
        vector<int> v(n);
        for(int j = 0; j < n; j++)
            fe >> v[j];
        d.push_back(v);
    }

    fe.close();
}


int main(int argc, char** argv){

    struct timespec cgt1,cgt2;
    vector<int> solution;
    double ncgt;
    int max_generations, result;
    double score;
    string algorithm_name;

    if (argc < 3) {
        cout << "[ERROR] USO: <" << argv[0] << "> <fichero> <opcion>\n";
        cout << "Para elegir una opción seleccione un número:" << endl;
        cout << " 1: Algoritmo genético estándar." << endl;
        cout << " 2: Algoritmo genético baldwiniano." << endl;
        cout << " 3: Algoritmo genético lamarckiano." << endl;
        exit(-1);
    }

    // Read flow and distances
    input(argv[1]);
    // Get option
    string option = argv[2];

    // Get num of generations
    if(argc == 4){
        max_generations = atoi(argv[3]);
    }
    else{
        max_generations = 100;
    }

    // Start clock
    clock_gettime(CLOCK_REALTIME,&cgt1);

    // Execute algorithm

    switch(option[0]){
        case '1':
            algorithm_name = "ESTANDAR";
            solution = Generational(algorithm_name, max_generations);
        break;
        case '2':
            algorithm_name = "BALDWINIANO";
            solution = Generational(algorithm_name, max_generations);
        break;
        case '3':
            algorithm_name = "LAMARCKIANO";
            solution = Generational(algorithm_name, max_generations);
        break;
        case '4':
            algorithm_name = "ESTACIONARIO";
            solution = solution = Stationary(max_generations);
        break;
    }

    // Stop clock
    clock_gettime(CLOCK_REALTIME,&cgt2);

    // Get results
    ncgt = (double) (cgt2.tv_sec-cgt1.tv_sec)+(double) ((cgt2.tv_nsec-cgt1.tv_nsec)/(1.e+9));
    result = costFunction(solution);
    score = max(5.0-100.0*(costFunction(solution)-44759294.0)/44759294.0, 0.0);

    // Print results
    cout << algorithm_name << ": " << result << endl;
    cout << "Time:   " << ncgt << endl << endl;
    cout << "Nota: " << score << endl;

    for(int i = 0; i < n; i++)
        cout << solution[i] << " ";
    cout << endl << endl;

}
