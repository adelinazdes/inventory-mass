﻿
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
    double massa_section; //масса dm

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


    double cross_sectional_area() const {
        double D = get_inner_diameter();
        double S = M_PI * D * D / 4;
        return S;
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

//  структура для Проблемно-ориентированного буфера
struct density_viscosity_pressure_layer
{
    vector<double> density;
    vector<double> viscosity;
    vector<double> pressure;
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



vector<double>  euler(pipe myPipe, ring_buffer_t<density_viscosity_pressure_layer>& buffer, int i) {
    double p_0 = myPipe.p_0;
    double p_begin = myPipe.p_0;
    density_viscosity_pressure_layer& current_layer = buffer.current();
    density_viscosity_pressure_layer& previous_layer = buffer.previous();

    //вектор хранит значения давлений на первом слое
    vector<double> pressure_first_layer;
    if (i == 0) {
        size_t b = current_layer.density.size();
        for (size_t j = 0; j < b; j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / previous_layer.viscosity[j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer.pressure[j] = p_0;
            double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer.density[j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer.density[j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
            pressure_first_layer.push_back(previous_layer.pressure[j]); //добавление в вектор давления на первом слое
        }
    }

    if (i != 0) {
        for (size_t j = 0; j < current_layer.density.size(); j++) {
            double Re = myPipe.V * myPipe.get_inner_diameter() / previous_layer.viscosity[j];
            double lambda = hydraulic_resistance_isaev(Re, myPipe.get_relative_roughness());
            previous_layer.pressure[j] = p_0;
            double p_rachet = p_0 - myPipe.get_dx() * (lambda / myPipe.get_inner_diameter() * previous_layer.density[j] * pow(myPipe.V, 2) / 2 - M_G * previous_layer.density[j] * (myPipe.z_L - myPipe.z_0) / ((myPipe.N - 1) * myPipe.get_dx()));
            p_0 = p_rachet;
        }
    }

    return pressure_first_layer;
}

// запись в файл excel
void excel(const string& filename, pipe myPipe, density_viscosity_pressure_layer& previous_layer, int i, vector <double>& data, vector <double>& massa, vector <double>& wall, vector <double>& density) {

    for (int j = 0; j < myPipe.N; j++) {
        if (i == 0 && j == 0) {
            ofstream outFile(filename, ios::out);
            setlocale(LC_ALL, "ru");
            outFile << "время,координата,плотность, вязкость, давление, разность давления, масса, D, плотность рабочая " << "\n";
            outFile << myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer.density[j] << "," << previous_layer.viscosity[j] << "," << previous_layer.pressure[j] << "," << previous_layer.pressure[j] - data[j] << "," << massa[j]<< "," <<pow(4*wall[j]/ M_PI,0.5) << "," << density[j] << "\n";
            outFile.close();
        }
        else {
            ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
            outFile << myPipe.time << "," << (j)*myPipe.get_dx() << "," << previous_layer.density[j] << "," << previous_layer.viscosity[j] << "," << previous_layer.pressure[j] << "," << previous_layer.pressure[j] - data[j] << "," << massa[j] << "," << pow(4 * wall[j] / M_PI, 0.5) << "," << density[j] << "\n";
            outFile.close();

        }

    }
}

/// запись в файл excel значений давлений в последней точке каждого слоя
void excel_last_pressure(const string& filename, pipe myPipe, density_viscosity_pressure_layer& previous_layer, int i, vector <double>& wall_layer) {

    //давление в последней точке трубопровода на каждом слое
    if (i == 0) {
        ofstream outFile(filename, ios::out);
        setlocale(LC_ALL, "ru");
        outFile << "нули,время,давление, объем" << "\n";
        double volume=0;
        for (int j = 1; j < myPipe.N; j++) {
            double volume_section = wall_layer[j] * myPipe.get_dx();
            volume = volume + volume_section;
        }
        
        outFile << 0 << "," << myPipe.time << "," << previous_layer.pressure[myPipe.N - 1] << "," << volume << "\n";
        outFile.close();
    }
    else {
        double volume = 0;
        for (int j = 1; j < myPipe.N; j++) {
            double volume_section = wall_layer[j] * myPipe.get_dx();
            volume = volume + volume_section;
        }

        ofstream outFile(filename, ios::app); // Используйте ios::app для добавления данных
        outFile << 0 << "," << myPipe.time << "," << previous_layer.pressure[myPipe.N - 1] << "," << volume << "\n";
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

// @brief Ступенчатая партия, вытеснение
//TEST(BLOCK3, TASK1)
//{
//    pipe myPipe = block3_pipe();
//
//    vector<double> ro;
//    ro = { 800 }; //значения плотности на входе трубы
//
//    vector<double> u;
//    u = { 10e-6 }; //значения кинемат. вязкости на входе трубы
//
//    vector<double> pressure;
//    pressure = { 0 };//значения давлений
//
//    Проблемно-ориентированный буфер из 3 слоев
//    density_viscosity_pressure_layer initialization_layer;
//    initialization_layer.density = vector<double>(myPipe.N,900);
//    initialization_layer.viscosity = vector<double>(myPipe.N, 15e-6);
//    initialization_layer.pressure = vector<double>(myPipe.N,0);
//
//    ring_buffer_t<density_viscosity_pressure_layer> buffer(2, initialization_layer);
//
//    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 
//
//    for (size_t i = 0; i < myPipe.get_n() + 2; i++) {
//        if (i == 0) {
//            myPipe.time = myPipe.get_dt() * i;
//            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
//            party_layer(myPipe, ro[0], buffer.current().density, buffer.previous().density);
//            party_layer(myPipe, u[0], buffer.current().viscosity, buffer.previous().viscosity);
//            pressure_first_layer = euler(myPipe, buffer, i);
//            excel("3block_1.1.csv", myPipe, buffer.previous(), i, pressure_first_layer);
//            buffer.advance(1);
//        }
//        else {
//            myPipe.time = myPipe.get_dt() * i;
//            euler(myPipe, buffer, i);
//            excel_last_pressure("pressure_last_1.1.csv", myPipe, buffer.current(), i);
//            party_layer(myPipe, ro[0], buffer.current().density, buffer.previous().density);
//            party_layer(myPipe, u[0], buffer.current().viscosity, buffer.previous().viscosity);
//            excel("3block_1.1.csv", myPipe, buffer.previous(), i, pressure_first_layer);
//            buffer.advance(1);
//        }
//    }
//
//
//}


///// @brief Импульсная партия
//TEST(BLOCK3, TASK2) {
//
//    pipe myPipe = block3_pipe();
//
//    double impulse_party_duration = 5 * 3600; // длительность захода импульсной партии, сек
//    size_t impulse_party_layers = ceil(impulse_party_duration / myPipe.get_dt());// кол-во слоев для импульсной партии
//
//    // время моделирования - время вытеснения начальной партии плюс длительность импульсной партии
//    size_t total_layers = myPipe.get_n() + impulse_party_layers + 2;
//
//    vector<double> ro;
//    ro = {};
//    // заполняем значениями импульсной партии
//    for (int k = 0; k < impulse_party_layers; ++k) {
//        ro.push_back(870);
//    }
//
//    // Заполняем оставшиеся 50 элементов значением основной партии
//    for (int l = impulse_party_layers; l < total_layers; ++l) {
//        ro.push_back(830);
//    }
//
//    vector<double> u;
//    u = {  };
//    // заполняем значениями импульсной партии
//    for (int k = 0; k < impulse_party_layers; ++k) {
//        u.push_back(19e-6);
//    }
//
//    // Заполняем оставшиеся 50 элементов значением основной партии
//    for (int l = impulse_party_layers; l < total_layers; ++l) {
//        u.push_back(15e-6);
//    }
//    //Проблемно-ориентированный буфер из 3 слоев
//    density_viscosity_pressure_layer initialization_layer;
//    initialization_layer.density = vector<double>(myPipe.N, 830);
//    initialization_layer.viscosity = vector<double>(myPipe.N, 15e-6);
//    initialization_layer.pressure = vector<double>(myPipe.N, 0);
//
//    ring_buffer_t<density_viscosity_pressure_layer> buffer(2, initialization_layer);
//    
//    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 
//
//    for (size_t i = 0; i < total_layers; i++) {
//        if (i == 0) {
//            myPipe.time = myPipe.get_dt() * i;
//            euler(myPipe, buffer, i);
//            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.current(), i);
//            party_layer(myPipe, ro[i], buffer.current().density, buffer.previous().density);
//            party_layer(myPipe, u[i], buffer.current().viscosity, buffer.previous().viscosity);
//            pressure_first_layer = euler(myPipe, buffer, i);
//            excel("3block_1.2.csv", myPipe, buffer.previous(), i, pressure_first_layer);
//            buffer.advance(1);
//        }
//        else {
//            myPipe.time = myPipe.get_dt() * i;
//            euler(myPipe, buffer, i);
//            excel_last_pressure("pressure_last_1.2.csv", myPipe, buffer.current(), i);
//            party_layer(myPipe, ro[i], buffer.current().density, buffer.previous().density);
//            party_layer(myPipe, u[i], buffer.current().viscosity, buffer.previous().viscosity);
//            excel("3block_1.2.csv", myPipe, buffer.previous(), i, pressure_first_layer);
//            buffer.advance(1);
//        }
//    }
//}


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



vector<double> time_step_characteristics_method(pipe myPipe, double& input_Q, int i, vector<double> time_modeling) {

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
    double step = time_modeling / step_diskretizatsii;
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




// плотность с учетом сжимаемости

vector<double> density(pipe myPipe, density_viscosity_pressure_layer& previous_layer) {

    vector<double> density(myPipe.N, 0.0);
    double dx = myPipe.get_dx();
    double p_0 = 101325;
    double koeff_density = 1.4*pow(10, 9);
    for (int j = 1; j < myPipe.N; j++) {
        double density_step = previous_layer.density[j] * (1 +(previous_layer.pressure[j] - p_0)/ koeff_density);
        density[j] = density_step;
    }
    return density;
}




//расчет массы в каждой координатной точке  сетки

vector<double> massa (pipe myPipe,  density_viscosity_pressure_layer& previous_layer, vector<double> wall_layer, vector<double> density_layer) {
    
    vector<double> massa(myPipe.N, 0.0);
    double dx = myPipe.get_dx();

    for (int j = 1; j < myPipe.N; j++) {
        double mean_density = 0.5 * (density_layer[j - 1] + density_layer[j]);
        double mean_S = 0.5 * (wall_layer[j - 1] + wall_layer[j]);
        double massa_section = dx * mean_S * mean_density;
        massa[j] = massa[j - 1] + massa_section;
    }
    return massa;
}



//расчет S в каждой координатной точке  сетки с учетом растяжимости стенки трубопровода 

vector<double> wall (pipe myPipe, density_viscosity_pressure_layer& previous_layer) {

    vector<double> wall (myPipe.N, 0.0);
    double S = myPipe.cross_sectional_area();
    double dx = myPipe.get_dx();
    double E = (2.1) * pow(10, 11);
    double p_0 = 101325;
    double v = 0.28;
    double koeff_wall = S * (1 - pow(v, 2)) / E / myPipe.b;
    for (int j = 1; j < myPipe.N; j++) {
        double wall_step = S * (1+koeff_wall*(previous_layer.pressure[j] - p_0));
        wall[j] = wall_step;
    }
    return wall;
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

    //вектор состоящий из зна-ий времени моделирования (шаг метода характеристик)
    vector<double> time_modeling;
    time_modeling.push_back(0.0);

    vector<double> Q_int;
    Q_int.push_back(Q[0]);

    for (size_t i = 0; i < total_layers - 1; i++) {
        time_modeling = time_step_characteristics_method(myPipe, Q_int[i], i, time_modeling);
        Q_int.push_back(interpolation_scalar(time_modeling[i + 1], Q));

    }

    //вектор состоящий из зна-ий при интепроляции
    vector<double> ro_int = interpolation_vector(time_modeling, ro);
    vector<double> u_int = interpolation_vector(time_modeling, u);
    vector<double> pressure_int = interpolation_vector(time_modeling, pressure);
    
    vector<double> pressure_first_layer;// вектор содержащий значения давлений первого слоя 

    //Проблемно-ориентированный буфер из 3 слоев
    density_viscosity_pressure_layer initialization_layer;
    initialization_layer.density = vector<double>(myPipe.N, ro[0]);
    initialization_layer.viscosity = vector<double>(myPipe.N, u[0]);
    initialization_layer.pressure = vector<double>(myPipe.N, pressure[0]);

    ring_buffer_t<density_viscosity_pressure_layer> buffer(2, initialization_layer);

    for (size_t i = 0; i < total_layers; i++) {
        myPipe.Q = Q_int[i];
        myPipe.V = myPipe.get_V();
        myPipe.p_0 = pressure_int[i];
        myPipe.time = time_modeling[i];
        euler(myPipe, buffer, i);
        //excel_last_pressure("pressure_last_1.3.csv", myPipe, buffer.previous(), i, wall_layer);
        party_layer(myPipe, ro_int[i], buffer.current().density, buffer.previous().density);
        party_layer(myPipe, u_int[i], buffer.current().viscosity, buffer.previous().viscosity);
        
        if (i == 0) {
            pressure_first_layer = euler(myPipe, buffer, i);
        }
        
        vector <double> wall_layer = wall(myPipe, buffer.previous());
        vector <double> density_layer = density(myPipe, buffer.previous());
        vector <double> massa_layer =massa(myPipe, buffer.previous(), wall_layer,density_layer);
        excel_last_pressure("pressure_last_1.3.csv", myPipe, buffer.previous(), i, wall_layer);
        excel("3block_1.3.csv", myPipe, buffer.previous(), i, pressure_first_layer, massa_layer, wall_layer, density_layer);

        buffer.advance(1);
    }

}

