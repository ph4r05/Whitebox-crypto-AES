#ifndef NTL_matrix__H
#define NTL_matrix__H

#include <NTL/tools.h>
#include <NTL/vector.h>


// matrix templates

NTL_OPEN_NNS


template<class T> 
class Mat {  
public:  
  
   // pseudo-private fields
   Vec< Vec<T> > _mat__rep;  
   long _mat__numcols;  



   // really public fields

   typedef typename Vec<T>::value_type value_type;
   typedef typename Vec<T>::reference reference;
   typedef typename Vec<T>::const_reference const_reference;
  
  
   Mat() : _mat__numcols(0) { }  
   Mat(const Mat<T>& a);  
   Mat& operator=(const Mat<T>& a);  
   ~Mat() { }  
  
   Mat(INIT_SIZE_TYPE, long n, long m);  
  
   void kill();  
  
   void SetDims(long n, long m);  
  
   long NumRows() const { return _mat__rep.length(); }  
   long NumCols() const { return _mat__numcols; }  
  
   Vec<T>& operator[](long i) { return _mat__rep[i]; }  
   const Vec<T>& operator[](long i) const { return _mat__rep[i]; }  
  
   Vec<T>& operator()(long i) { return _mat__rep[i-1]; }  
   const Vec<T>& operator()(long i) const { return _mat__rep[i-1]; }  
  
   reference operator()(long i, long j) { return _mat__rep[i-1][j-1]; }  
   const_reference operator()(long i, long j) const   
      { return _mat__rep[i-1][j-1]; }  

   const_reference get(long i, long j) const { return _mat__rep[i].get(j); }
   void put(long i, long j, const T& a) { _mat__rep[i].put(j, a); }

   template <class U>
   void put(long i, long j, const U& a) { _mat__rep[i].put(j, a); }

  
   long position(const Vec<T>& a) const { return _mat__rep.position(a); } 
   long position1(const Vec<T>& a) const { return _mat__rep.position1(a); } 
   Mat(Mat<T>& x, INIT_TRANS_TYPE) :  
    _mat__rep(x._mat__rep, INIT_TRANS), _mat__numcols(x._mat__numcols) { }  
};  
 
template<class T> 
inline const Vec< Vec<T> >& rep(const Mat<T>& a)  
   { return a._mat__rep; }  
  

template<class T>
Mat<T>::Mat(const Mat<T>& a) : _mat__numcols(0)  
{  
   SetDims(a.NumRows(), a.NumCols());  
   _mat__rep = a._mat__rep;  
}  
  
template<class T>
Mat<T>& Mat<T>::operator=(const Mat<T>& a)  
{  
   SetDims(a.NumRows(), a.NumCols());  
   _mat__rep = a._mat__rep;  
   return *this;  
}  
  
template<class T>
Mat<T>::Mat(INIT_SIZE_TYPE, long n, long m) : _mat__numcols(0)
{  
   SetDims(n, m);  
}  
  
template<class T>
void Mat<T>::kill()  
{  
   _mat__numcols = 0;  
   _mat__rep.kill();  
}  
  
template<class T>
void Mat<T>::SetDims(long n, long m)  
{  
   if (n < 0 || m < 0)  
      Error("SetDims: bad args");  
  
   if (m != _mat__numcols) {  
      _mat__rep.kill();  
      _mat__numcols = m;  
   }  
        
   long oldmax = _mat__rep.MaxLength();  
   long i;  
   _mat__rep.SetLength(n);  
  
   for (i = oldmax; i < n; i++)  
      _mat__rep[i].FixLength(m);  
}  
     
        
template<class T>
void MakeMatrix(Mat<T>& x, const Vec< Vec<T> >& a)  
{  
   long n = a.length();  
  
   if (n == 0) {  
      x.SetDims(0, 0);  
      return;  
   }  
  
   long m = a[0].length();  
   long i;  
  
   for (i = 1; i < n; i++)  
      if (a[i].length() != m)  
         Error("nonrectangular matrix");  
  
   x.SetDims(n, m);  
   for (i = 0; i < n; i++)  
      x[i] = a[i];  
}  
  
template<class T>
void swap(Mat<T>& X, Mat<T>& Y)  
{  
   swap(X._mat__numcols, Y._mat__numcols);  
   swap(X._mat__rep, Y._mat__rep);  
}  
  
template<class T>
long operator==(const Mat<T>& a, const Mat<T>& b)  
{  
   if (a.NumCols() != b.NumCols())  
      return 0;  
  
   if (a.NumRows() != b.NumRows())  
      return 0;  
  
   long n = a.NumRows();  
   long i;  
  
   for (i = 0; i < n; i++)  
      if (a[i] != b[i])  
         return 0;  
  
   return 1;  
}  
  
  
template<class T>
long operator!=(const Mat<T>& a, const Mat<T>& b)  
{  
   return !(a == b);  
}  


template<class T>
NTL_SNS istream& operator>>(NTL_SNS istream& s, Mat<T>& x)  
{  
   Vec< Vec<T> > buf;  
   s >> buf;  
   MakeMatrix(x, buf);  
   return s;  
}  
  
template<class T>
NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const Mat<T>& a)  
{  
   long n = a.NumRows();  
   long i;  
   s << "[";  
   for (i = 0; i < n; i++) {  
      s << a[i]; 
      s << "\n"; 
   }  
   s << "]";  
   return s;  
}  


// conversion

template<class T, class S>
void conv(Mat<T>& x, const Mat<S>& a)
{  
   x.SetDims(a.NumRows(), a.NumCols());  
   conv(x._mat__rep, a._mat__rep);  
}  



NTL_CLOSE_NNS


#endif
