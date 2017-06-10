//#include "Vector3.h"
//
//template<typename T>
//
//Vector3::Vector3 : x(T(0)), y(T(0)), z(T(0))()
//{
//}
//
//Vector3::Vector3(T xx) : x(xx), y(xx), z(xx)()
//{
//}
//
//Vector3::Vector3(T xx, T yy, T zz) : x(xx), y(yy), z(zz)
//{
//}
//
//Vector3::~Vector3()
//{
//}
//
//Vector3 & Vector3<T>::normalize() //& antes do método faz com que ele mesmo seja usado como parametro. O nome disto é extension methods.
//{
//	T nor2 = length2();
//	if (nor2 > 0) {
//		T invNor = 1 / sqrt(nor2);
//		x *= invNor, y *= invNor, z *= invNor;
//	}
//	return *this;
//}
//
//Vector3<T> Vector3<T>::operator*(const T & p_f) const //Lembrete: As sugestões do visual studio no .h para delcarar métodos no cpp são PERFEITAS, poupam tempo e bugs.
//{
//	return Vector3<T>(x * p_f, y * p_f, z * p_f);
//}
//
//Vector3<T> Vector3<T>::operator*(const Vector3<T>& p_v) const
//{
//	return Vector3<T>(x * p_v.x, y * p_v.y, z * p_v.z);
//}
//
//T Vector3<T>::dot(const Vector3<T>& v) const
//{
//	return x * v.x + y * v.y + z * v.z;
//}
//
//Vector3<T> Vector3<T>::operator-(const Vector3<T>& p_v) const
//{
//	return Vector3<T>(x - p_v.x, y - p_v.y, z - p_v.z);
//}
//
//Vector3<T> Vector3<T>::operator+(const Vector3<T>& p_v) const
//{
//	return Vector3<T>(x + p_v.x, y + p_v.y, z + p_v.z);
//}
//
//Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& p_v)
//{
//	x += p_v.x, y += p_v.y, z += p_v.z; return *this;
//}
//
//Vector3<T>& Vector3<T>::operator*=(const Vector3<T>& p_v)
//{
//	x *= p_v.x, y *= p_v.y, z *= p_v.z; return *this;
//}
//
//Vector3<T> Vector3<T>::operator-() const
//{
//	return Vector3<T>(-x, -y, -z);
//}
//
//T Vector3<T>::length2() const
//{
//	return x * x + y * y + z * z;
//}
//
//T Vector3<T>::length() const
//{
//	return sqrt(length2());
//}
//
//friend ostream & operator << (ostream &p_os, const Vector3<T> &p_v)
//{
//	p_os << "[" << p_v.x << " " << p_v.y << " " << p_v.z << "]";
//	return p_os;
//}