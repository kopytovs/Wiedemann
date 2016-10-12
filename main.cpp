#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <list>
#include <ctime>


using namespace std;

ifstream in("/Users/Kopytov/Documents/AppleScript/Wideman_alg/Wideman_alg/input.txt");
ofstream on("/Users/Kopytov/Documents/AppleScript/Wideman_alg/Wideman_alg/output.txt");

#define cin in
#define cout on

//int n;
//vector < string > simple_rehash;
//vector <vector< string >> set_hash;

//int cnt_str;
//int cnt_requests;


//Алгоритм Берлекэмпа-Месси
list<int> messi(vector<int> mas);

//Перемножение матриц
vector< vector<int> > mult_matrix(vector< vector<int> > m1, vector< vector<int> > m2 );

// Суммирование матриц
vector< vector<int> > sum_matrix(vector< vector<int> > m1, vector< vector<int> > m2 );

// Есть ли хоть один не нулевой элемент в матрице
bool is_zero(vector< vector<int> > m1 );

vector< vector <int> > A;
vector< vector <int> > b;
vector< vector <int> > bx;
vector< vector <int> > y;
vector< vector <int> > u;

vector< vector <int> > A0;

vector< vector< vector <int> > > exp_of_A;

int main()
{
    int n;
    cin>>n;
    A.resize(n);
    A0.resize(n);
    b.resize(n);
    y.resize(n);
    exp_of_A.clear();
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            int t;
            cin >> t;
            A[i].push_back(t);
            if( i == j )
                A0[i].push_back(1);
            else
                A0[i].push_back(0);
            
        }
    }
    for(int i= 0 ; i< n; i++)
    {
        int t;
        cin >> t;
        b[i].push_back(t);
        y[i].push_back(0);
    }
    
    exp_of_A.push_back(A0);
    
    // Возводим матрицу в степень
    for(int i= 1; i <= 2*n; i++)
    {
        exp_of_A.push_back( mult_matrix(A,exp_of_A[i-1]) );
    }
    
    unsigned long q = clock();
    
    // Алгоритм Видемана
    bx = b;
    int d = 0;
    
    srand((unsigned)time(0));
    
    // Проверяем, не нашли ли мы решение. Для этого проверяем измененённую правую часть (1 шаг)
    while (!is_zero(bx))
    {
        u.clear();
        // формируем случайный вектор (2 шаг)
        while (is_zero(u))
        {
            u.clear();
            u.resize(1);
            for(int i =0; i< n ; i++ )
            {
                
                u[0].push_back(rand() % 2);
                //u[0].push_back(1);
            }
        }
        
        vector<int> z;
        list<int> f;
        z.clear();
        f.clear();
        // В этом случае мы зацикливаемся, так как не можем получить последовательность на 4 шаге
        // А это может быть, только когда нет решений
        if (d >= n)
        {
            cout << "Определитель матрицы нулевой!" << endl;
            break;
        }
        
        
        // Вычисляем первые 2(n-d) последовательности( 4 шаг )
        for(int i= 0; i < 2*(n-d); i++ ) 
        { 
            vector< vector <int> > t; 
            t = mult_matrix( exp_of_A[i],bx); 
            t = mult_matrix( u , t ); 
            z.push_back(t[0][0]); 
        } 
        // Вычисляем минимальный многочлен последовательности(5 шаг) 
        f = messi(z); 
        // нормируем минимальный многочлен 
        f.pop_back(); 
        
        int step; 
        vector< vector <int> > t; 
        t.clear(); 
        t.resize(n); 
        for(int j=0; j<n; j++) 
            t[j].resize(n); 
        step=0; 
        //for (std::list<int>::iterator it=f.begin(); it != f.end(); ++it) 
        // Степени матрицы лежат задом на перед 
        for (std::list<int>::reverse_iterator it=f.rbegin(); it!=f.rend(); ++it) 
        { 
            if (*it == 1 ) 
                t = sum_matrix(t,exp_of_A[step]); 
            step++; 
        } 
        /*if (step > 0) 
         step--;*/ 
        
        // Переприсваеваем переменные(шаг 6) 
        t = mult_matrix(t,bx); 
        y= sum_matrix(y,t); 
        t.clear(); 
        t.resize(n); 
        for(int j=0; j<n; j++) 
            t[j].resize(n); 
        
        t = mult_matrix(A,y); 
        bx = sum_matrix(b,t); 
        d += step; 
    } 
    
    
    
    
    // Выводим ответ 
    for (int i = 0; i < y.size(); i++) 
    { 
        for(int j = 0; j < y[i].size(); j++) 
        { 
            cout << " " << y[i][j];
        } 
        cout<<endl;
    } 
    q = (clock()-q);
    cout << q << "msec" << endl;
    return 0;
    
}

//Алгоритм Берлекэмпа-Месси
list<int> messi(vector<int> mas)
{
    
    list<int> res;// результирующая последовательность
    res.clear();
    int delta = 0;// невязка
    list<int> c;
    c.clear();
    c.push_back(1);// коэффициенты порождающей последовательности
    list<int> b;
    b.clear();
    b.push_back(1);// временная переменная, которая хранит модифицирующий вектор
    int m = 0;
    
    for(int i=0; i<mas.size() ; i++)
    {
        int tek_c_ind;
        tek_c_ind = 0;
        delta = 0;
        // находим невязку
        for (std::list<int>::iterator it=c.begin(); it != c.end(); ++it)
        {
            delta += mas[i-tek_c_ind] * *it;
            tek_c_ind++;
        }
        delta%=2;
        // Если нет невязки просто увеличиваем степень многочлена
        if (delta == 0)
        {
            b.push_front(0);
            m = i;
        }
        else
        {
            // Если невязка не нулевая
            list<int> t;// переменная в которой будет лежать x * b(x) + c(x)
            list<int> tmp_b; // просто временная переменная, нужна т.к. удаляем элементы при суммировании
            t.clear();
            tmp_b = b;
            //for( int j = 0; j < i - m; j++)
            tmp_b.push_front(0);
            list<int> tmp_c;// просто временная переменная, нужна т.к. удаляем элементы при суммировании
            tmp_c = c;
            //суммируем t и с
            while (!tmp_b.empty() || !tmp_c.empty())
            {
                int tek_sum;
                tek_sum = 0;
                if (!tmp_b.empty()) {
                    tek_sum += tmp_b.front();
                    tmp_b.pop_front();
                }
                if (!tmp_c.empty()) {
                    tek_sum += tmp_c.front();
                    tmp_c.pop_front();
                }
                tek_sum%=2;
                t.push_back(tek_sum);
            }
            list<int> tt;
            tt = t;
            // КОНЕЧНЫЕ НУЛИ ЗНАЧАЩИЕ ДЛЯ МНОГОЧЛЕНА!!!!! Их нельзя удалять!!
            /*   for (std::list<int>::reverse_iterator rit=tt.rbegin(); rit!=tt.rend(); ++rit)
             {
             if (*rit == 1) break;
             if (*rit == 0) t.pop_back();
             }*/
            
            b.push_front(0);
            // Этого нет в описанном алгоритме, поэтому описанный алгоритм в лабе не выводит минимальную
            // порождающую последовательность, а выводит любую порождающую последовательность.
            // Чтобы была минимальная надо на каждом шаге сравнивать длину
            // b*b(x) и x * b(x) + c(x). Что меньше - то и сохранять.
            if (c.size() <= b.size())
                b = c;
            c = t;
        }
    }
    res = c;
    return res;
}

//Перемножение матриц
vector< vector<int> > mult_matrix(vector< vector<int> > m1, vector< vector<int> > m2 )
{
    vector< vector<int> > res;
    res.clear();
    
    unsigned long int max_j = m2[0].size();
    unsigned long int max_k = m2.size();
    res.resize(m1.size());
    
    for( int i =0; i < m1.size(); i++)
    {
        for (int j = 0; j < max_j; j++ )
        {
            int ssum = 0;
            
            for (int k = 0; k < max_k; k++)
                ssum += m1[i][k] * m2[k][j];
            
            res[i].push_back(ssum%2);
        }
    }
    return res;
}

// Суммирование матриц
vector< vector<int> > sum_matrix(vector< vector<int> > m1, vector< vector<int> > m2 )
{
    vector< vector<int> > res;
    res.clear();
    res.resize(m1.size());
    for(int i = 0; i < m1.size(); i++)
        for(int j = 0; j< m2[i].size(); j++)
            res[i].push_back( (m1[i][j] + m2[i][j]) %2 );
    
    return res;
}

// Есть ли хоть один не нулевой элемент в матрице
bool is_zero(vector< vector<int> > m1 )
{
    
    for(int i = 0; i < m1.size(); i++)
        for(int j = 0; j< m1[i].size(); j++)
            if (m1[i][j] > 0) return false;
    
    return true;
    
}

