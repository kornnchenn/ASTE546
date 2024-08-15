
/*
int a = 1;
void my_fun(int a) {
 a = -20;
}
int main() {
 int b = 0;
 b += a; // step 1
 int a = 2;
 for (int a, i = 0; i<10; i++) {
 a+=1;
 }
 b += a; // step 2
 my_fun(a);
 b += a; // step 3
}
*/
/*
void f1(int a) { a = a + 1; }
void f2(int* a) {*a = *a + 1; }
void f3(int &a) {a = a + 1;}
int main() {
  int a = 0;
  f1(a);
  f2(&a);
  f3(a);
}
*/

struct MyData {
  int x;
};
void f1(MyData m1, MyData m2) { m1.x = m2.x; }
void f2(MyData m1, MyData &m2) { m1.x = m2.x; }
void f3(MyData &m1, MyData m2) { m1.x = m2.x; }
void f4(MyData &m1, MyData &m2) { m1.x = m2.x; }
void f5(MyData &m1, const MyData &m2) { m1.x = m2.x; }
void f6(const MyData &m1, const MyData &m2) { m1.x = m2.x; }
int main() {
  MyData A, B;
  B.x = 10;
  f1(A,B);  
  f2(A,B), 
  f3(A,B);
  f4(A,B);  
  f5(A,B), 
  f6(A,B);
}