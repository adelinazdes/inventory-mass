#include <iomanip>
#include <iostream>
#include <clocale>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <locale.h>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>
#include <cmath>
#include <gtest/gtest.h>   

using namespace std;
using namespace pde_solvers;

struct pipe
{
    double L; // протяженность участка
    double V; //скорость нефти
    double D_vnesh; //диаметр внешний
    double b; //толщина стенки
    double abc; //абсолютная шероховатость в м
    double z_0, z_L; //высотные отметки начала,конца
    double ro; //плотность
    double u; // кинематическая вязкость Стоксы в м2/с
    double Q; //расход м3/c
    double p_0;// давление в начале участка
    double p_L;// давление в начале участка
    double lambda;//коэфф.гидравл.сопротивления
    double t_w;//касательное напряжение трения
    double h; //шаг по координате расчетной сетки, в метрах
    double Re;
    double N; //КОЛ-ВО ТОЧЕК


    double time;//шаг моделирования (по времени)

    double get_inner_diameter() const {
        return D_vnesh - 2 * b;
    }
    double get_relative_roughness() const {
        return abc / get_inner_diameter();
    }

    double get_inner_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return V * S;
    }

    double get_V() const {
        double D = get_inner_diameter();
        return 4 * Q / (3.1415 * pow(D, 2));
    }

    double get_Re() const {
        double D = get_inner_diameter();
        return get_V() * D / u;
    }
    //касательное напряжение трения
    double get_t_w() const {
        return lambda / 8 * ro * pow(get_V(), 2);
    }

    double get_dx() const {
        return L / (N - 1);
    }

    //шаг метода характеристик
    double get_dt() const {
        return get_dx() / V;
    }

    double get_n() const {
        return (L / V) / get_dt();
    }
};


// след. слой предыдущий слой + новое значение на входе
void party_layer(pipe myPipe, double parametr, vector<double>& current_layer, vector<double>& previous_layer) {

    //смещение предыдущего слоя и запись граничного условия
    current_layer[0] = parametr;
    for (size_t i = 1; i < myPipe.N; i++)
    {
        current_layer[i] = previous_layer[i - 1];
    }
}

vector<double>  euler(pipe myPipe, ring_buffer_t<vector<vector<double>>>& buffer, int i)  {
    double p_0 = myPipe.p_0;
    double p_begin = myPipe.p_0;
    vector<vector<double>>& current_layer = buffer.current();
    vector<vector<double>>& previous_layer = buffer.previous();
   
    //вектор хранит значения давлений на первом слое
    vector<double> pressure_first_layer;
    if (i == 0) {
        size_t b = current_layer[0].size();
        for (size_t j = 0; j < b; j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / previous_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
            double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            pressure_first_layer.push_back(previous_layer[2][j]); //добавление в вектор давления на первом слое
        }
    }

    if (i != 0) {
        for (size_t j = 0; j < current_layer[0].size(); j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / current_layer[1][j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer[2][j] = p_0;
             double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer[0][j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer[0][j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
        }
    }

    return pressure_first_layer;
}
// запись в файл excel
void excel(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i, vector <double>& data) {

    for (int j = 0; j < myPipe.N; j++) {
        if (i == 0 && j == 0) {
            ofstream outFile(filename, ios::out);
            setlocale(LC_ALL, "ru");
            outFile << "время,координата,плотность, вязкость, давление, разность давления" << "\n";
            outFile << myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();
        }
        else {
            ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
            outFile << myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer[0][j] << "," << previous_layer[1][j] << "," << previous_layer[2][j] << "," << previous_layer[2][j] - data[j] << "\n";
            outFile.close();

        }

    }
}
// запись в файл excel занчений давлений в последней точке каждого слоя
void excel_last_pressure(const string& filename, pipe myPipe, vector<vector<double>>& previous_layer, int i) {

//давление в последней точке трубопровода на каждом слое
    if (i == 0) {
        ofstream outFile(filename, ios::out);
        setlocale(LC_ALL, "ru");
        outFile << "нули,время,давление" << "\n";
        outFile << 0 << "," << myPipe.time << "," << previous_layer[2][myPipe.N - 1] << "\n";
        outFile.close();
    }
    else {
        ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
        outFile << 0 << "," << myPipe.time << "," << previous_layer[2][myPipe.N - 1] << "\n";
        outFile.close();

    }


}

/// @brief Труба из заданий блока 3
pipe block3_pipe()
{
    pipe myPipe;
    myPipe.p_0 = 6e6;
    myPipe.L = 100e3;
    myPipe.D_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_L = 50;
    myPipe.V = 0.5;
    myPipe.N = 101;
    myPipe.abc = 15e-6;
    return myPipe;
}

/// @brief Ступенчатая партия, вытеснение
TEST(BLOCK3, TASK1) 
{
    pipe myPipe = block3_pipe();
   

    vector<double> ro;
    ro = { 800 }; //значения плотности на входе трубы

    vector<double> u;
    u = { 10e-6 }; //значения кинемат. вязкости на входе трубы

    vector<double> pressure;
    pressure = { 0 };//значения давлений


    vector<double> ro_begin(myPipe.N, 900);
    vector<double> u_begin(myPipe.N, 15e-6);
    vector<double> pressure_begin(myPipe.N, 0);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });

    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 

   for (size_t i = 0; i < myPipe.get_n() + 2; i++) {
        if (i == 0) {
            myPipe.time = myPipe.get_dt()*i;
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
            party_layer(myPipe, ro[0], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[0], buffer.current()[1], buffer.previous()[1]);
            
            pressure_first_layer = euler(myPipe, buffer, i);
            excel("3block_1.1.csv",myPipe, buffer.previous(), i, pressure_first_layer);
            buffer.advance(1);
        }
        else {
            myPipe.time = myPipe.get_dt() * i;
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
            party_layer(myPipe, ro[0], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[0], buffer.current()[1], buffer.previous()[1]);
            excel("3block_1.1.csv",myPipe, buffer.previous(), i, pressure_first_layer);
            buffer.advance(1);
        }
    }


}


/// @brief Импульсная партия
TEST(BLOCK3, TASK2) {

    pipe myPipe = block3_pipe();
    
    double impulse_party_duration = 5*3600; // длительность захода импульсной партии, сек
    size_t impulse_party_layers = ceil(impulse_party_duration / myPipe.get_dt());// кол-во слоев для импульсной партии

    // время моделирования - время вытеснения начальной партии плюс длительность импульсной партии
    size_t total_layers = myPipe.get_n() + impulse_party_layers + 2; 

    vector<double> ro;
    ro = {};
    // заполняем значениями импульсной партии
    for (int k = 0; k < impulse_party_layers; ++k) {
        ro.push_back(990);
    }

    // Заполняем оставшиеся 50 элементов значением основной партии
    for (int l = impulse_party_layers; l < total_layers; ++l) {
        ro.push_back(900);
    }

    vector<double> u;
    u = {  };
    // заполняем значениями импульсной партии
    for (int k = 0; k < impulse_party_layers; ++k) {
        u.push_back(19e-6);
    }

    // Заполняем оставшиеся 50 элементов значением основной партии
    for (int l = impulse_party_layers; l < total_layers; ++l) {
        u.push_back(15e-6);
    }

    vector<double> ro_begin(myPipe.N, 900);
    vector<double> u_begin(myPipe.N, 15e-6);
    vector<double> pressure_begin(myPipe.N, 0);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });

    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 

    for (size_t i = 0; i < total_layers; i++) {
        if (i == 0) {
            myPipe.time = myPipe.get_dt() * i;
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.previous(), i);
            party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
            pressure_first_layer = euler(myPipe, buffer, i);
            excel("3block_1.2.csv", myPipe, buffer.previous(), i, pressure_first_layer);
            buffer.advance(1);
        }
        else {
            myPipe.time = myPipe.get_dt() * i;
            euler(myPipe, buffer, i);
            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.previous(), i);
            party_layer(myPipe, ro[i], buffer.current()[0], buffer.previous()[0]);
            party_layer(myPipe, u[i], buffer.current()[1], buffer.previous()[1]);
            excel("3block_1.2.csv", myPipe, buffer.previous(), i, pressure_first_layer);
            buffer.advance(1);
        }
    }
}


/// @brief Труба из заданий блока 3_ задача 3
pipe block3_task3_pipe()
{
    pipe myPipe;
    myPipe.L = 10e4;
    myPipe.D_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_L = 50;
    myPipe.N = 101;
    myPipe.abc = 15e-6;
    return myPipe;
}



vector<double> time_step_characteristics_method (pipe myPipe, double& input_Q, int i, vector<double> time_modeling) {
     
    myPipe.Q = input_Q;
    myPipe.V = myPipe.get_V();
   //шаг метода характеристик (зависимость от расхода)
    double step_characteristics_method = myPipe.get_dt();
    if (i == 0) {
            time_modeling.push_back(step_characteristics_method);
    }
    else {
        double time = time_modeling.back() + step_characteristics_method;
        time_modeling.push_back(time);
     }
  
    return time_modeling;
}

double interpolation_scalar(double time_modeling, vector<double>& input) {
    double interpolation_value;
    double step_diskretizatsii = 2000;
    double step = time_modeling/ step_diskretizatsii;
    double first_slice = floor(step);
    double second_slice = ceil(step);
        if (first_slice >= input.size()) {
            first_slice = input.size() - 1;
        }
        if (second_slice >= input.size()) {
            second_slice = input.size() - 1;
        }

        double x = time_modeling;
        double x1 = first_slice * step_diskretizatsii; //время дискретизации имеющегося предыдущего значения 
        double x2 = second_slice * step_diskretizatsii; //время дискретизации имеющегося следующего значения 
        double y1 = input[first_slice];
        double y2 = input[second_slice];

        double y;
        if (abs(x - x1) < 1e-4) {
            y = y1;
        }
        else {
            y = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }
           interpolation_value = y;
   return interpolation_value;
}



vector<double> interpolation_vector(vector<double> time_modeling, vector<double>& input) {
    vector<double> interpolation_value;
    double step_diskretizatsii = 2000;
    for (size_t i = 0; i < time_modeling.size(); i++) {
        interpolation_value.push_back(interpolation_scalar(time_modeling[i], input));
        
    }
    return interpolation_value;
}





/// @brief задача 3
TEST(BLOCK3, TASK3) {

    pipe myPipe = block3_task3_pipe();
    
    //шаг дискретизации краевых условий (время через которое нач.условия изменяется)
    double step_diskretizatsii = 2000;
    
    //шаг метода характеристик (зависимость от расхода)
    double step_characteristics_method;

    //значения плотности на входе трубы
    vector<double> ro;
    ro = { 900,
        880,
        880,
        890,
        890,
        880,
        880,
        870 };

    //значения кинемат. вязкости на входе трубы
    vector<double> u;
    u = { 0.000015,
        0.000013,
        0.000013,
        0.000014,
        0.000014,
        0.000013,
        0.000013,
        0.000012 }; 
   
    //значения давления на входе трубы
    vector<double> pressure;
    pressure = { 6000000,
        5800000,
        5800000,
        5900000,
        5900000,
        5800000,
        5800000,
        5700000 };

    //значения расхода
    vector<double> Q;
    Q = { 0.192325,
        0.200000,
        0.210000,
        0.200000,
        0.180000,
        0.210000,
        0.210000,
        0.210000 };

    //кол-во слоев = кол-во значений на входе
    size_t  total_layers = Q.size();

    //векторы для буфера
    vector<double> ro_begin(myPipe.N, ro[0]);
    vector<double> u_begin(myPipe.N, u [0]);
    vector<double> pressure_begin(myPipe.N, pressure[0]);
    
    //вектор состоящий из зна-ий времени моделирования (шаг метода характеристик)
    vector<double> time_modeling;
    time_modeling.push_back(0.0);
    vector<double> Q_int;


    vector<double> ro_int;
    vector<double> u_int;
    vector<double> pressure_int;

    Q_int.push_back(Q[0]);

    for (size_t i = 0; i < total_layers; i++) {
        time_modeling = time_step_characteristics_method(myPipe, Q_int[i], i, time_modeling);
        Q_int.push_back(interpolation_scalar(time_modeling[i + 1], Q));
             
    }

    //вектор состоящий из зна-ий при интепроляции
    ro_int = interpolation_vector(time_modeling, ro);
    u_int = interpolation_vector(time_modeling, u);
    pressure_int = interpolation_vector(time_modeling, pressure);

    ring_buffer_t<vector<vector<double>>> buffer(2, { ro_begin, u_begin,pressure_begin });

    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 
        
    for (size_t i = 0; i < total_layers; i++) {
        myPipe.Q = Q_int[i];
        myPipe.V = myPipe.get_V();
        myPipe.p_0 = pressure_int[i];
        myPipe.time = time_modeling[i];
        euler(myPipe, buffer, i);
        excel_last_pressure("pressure_last_1.3.csv", myPipe, buffer.previous(), i);
        party_layer(myPipe, ro_int[i], buffer.current()[0], buffer.previous()[0]);
        party_layer(myPipe, u_int[i], buffer.current()[1], buffer.previous()[1]);
        if (i == 0) {
            pressure_first_layer = euler(myPipe, buffer, i);
        }
        excel("3block_1.3.csv", myPipe, buffer.previous(), i, pressure_first_layer);
        buffer.advance(1);
          
     }
        







}
