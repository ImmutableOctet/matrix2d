#pragma once

#ifndef IMMUTABLEOCTET_MATRIX2D_HPP

#define IMMUTABLEOCTET_MATRIX2D_HPP

#include <type_traits>
#include <utility>

#include <cstddef>
#include <cassert>
#include <cmath>

//#define IMMUTABLEOCTET_MATRIX2D_TYPE_PUNNING

namespace immutableoctet
{
	template <typename ValueType=float, typename VectorType=void>
	struct basic_matrix_2d
	{
		public:
			using value_type = ValueType;

			static constexpr bool has_vector_type = (!std::is_same<VectorType, void>::value); // inline

			static constexpr basic_matrix_2d identity()
			{
				return basic_matrix_2d {};
			}

			static constexpr std::size_t matrix_size()
			{
				return 6; // sizeof(raw_values);
			}

			constexpr basic_matrix_2d
			(
				value_type a, value_type b, value_type c, value_type d,
				value_type tx, value_type ty
			)
				: mA(a), mB(b), mC(c), mD(d), mTx(tx), mTy(ty)
			{
				//load_identity();
			}

			constexpr basic_matrix_2d() : basic_matrix_2d
			(
				static_cast<value_type>(1.0), static_cast<value_type>(0.0), static_cast<value_type>(0.0), static_cast<value_type>(1.0),
				static_cast<value_type>(0.0), static_cast<value_type>(0.0)
			)
			{
				//load_identity();
			}

			basic_matrix_2d(const basic_matrix_2d&) = default;
			basic_matrix_2d(basic_matrix_2d&&) noexcept = default;

			basic_matrix_2d& operator=(const basic_matrix_2d&) = default;
			basic_matrix_2d& operator=(basic_matrix_2d&&) noexcept = default;
			
#ifdef IMMUTABLEOCTET_MATRIX2D_TYPE_PUNNING
			constexpr value_type* data()
			{
				return raw_values;
			}

			constexpr const value_type* data() const
			{
				return raw_values;
			}
#endif

			constexpr std::size_t size() const
			{
				return matrix_size();
			}
		
			constexpr basic_matrix_2d concatenate(value_type mA, value_type mB, value_type mC, value_type mD, value_type mTx, value_type mTy) const
			{
				const auto a  = (mA * this->mA  + mC * this->mB);
				const auto b  = (mB * this->mA  + mD * this->mB);
				const auto c  = (mA * this->mC  + mC * this->mD);
				const auto d  = (mB * this->mC  + mD * this->mD);
				const auto tx = (mA * this->mTx + mC * this->mTy + mTx);
				const auto ty = (mB * this->mTx + mD * this->mTy + mTy);
			
				return basic_matrix_2d { a, b, c, d, tx, ty };
			}

			// Creates a matrix by concatenating this matrix with the input (`m`), combining their geometric effects.
			constexpr basic_matrix_2d concatenate(const basic_matrix_2d& m) const
			{
				return concatenate(m.mA, m.mB, m.mC, m.mD, m.mTx, m.mTy);
			}
		
			// Translates a copy of this matrix along the X and Y axes.
			constexpr basic_matrix_2d translate(value_type x, value_type y) const
			{
				auto m = *this;

				m.mTx += x;
				m.mTy += y;

				return m;
			}

			// Translates a copy of this matrix along the X axis.
			constexpr basic_matrix_2d translate_x(value_type x) const
			{
				auto m = *this;

				m.mTx += x;

				return m;
			}

			// Translates a copy of this matrix along the Y axis.
			constexpr basic_matrix_2d translate_y(value_type y) const
			{
				auto m = *this;

				m.mTy += y;

				return m;
			}
	
			// Applies a scaling transformation to a copy of this matrix.
			constexpr basic_matrix_2d scale(value_type sx, value_type sy) const
			{
				auto m = *this;

				m.mA *= sx;
				m.mB *= sy;
				m.mC *= sx;
				m.mD *= sy;
				m.mTx *= sx;
				m.mTy *= sy;

				return m;
			}
		
			// Applies a uniform scaling transformation to a copy of this matrix.
			constexpr basic_matrix_2d scale(value_type scalar) const
			{
				return scale(scalar, scalar);
			}
	
			// Applies a rotation on a copy of this matrix.
			// The angle specified is in radians.
			constexpr basic_matrix_2d rotate(value_type angle) const
			{
				return rotate(std::sin(angle), std::cos(angle));
			}
		
			// Applies a rotation on a copy of this matrix using sine and cosine values.
			constexpr basic_matrix_2d rotate(value_type r_sin, value_type r_cos) const
			{
				return rotate(r_sin, r_cos, -r_sin);
			}
		
			// Applies a rotation on a copy of this matrix using sine and cosine values.
			constexpr basic_matrix_2d rotate(value_type r_sin, value_type r_cos, value_type negative_r_sin) const
			{
				return concatenate(r_cos, r_sin, negative_r_sin, r_cos, 0.0, 0.0);
			}
		
			// Applies a shearing transformation to a copy of this matrix.
			constexpr basic_matrix_2d shear(value_type sx, value_type sy) const
			{
				auto m = *this;

				m.mC += sx;
				m.mB += sy;

				return m;
			}
	
			// Applies a uniform shearing transformation to a copy of this matrix.
			constexpr basic_matrix_2d shear(value_type scalar) const
			{
				auto m = *this;

				m.mC *= scalar; // +=
				m.mB *= scalar; // +=

				return m;
			}
	
			// Resets this matrix to an identity state.
			constexpr basic_matrix_2d& reset()
			{
				*this = basic_matrix_2d { 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 };

				return *this;
			}
		
			// Calculates the determinate value of the matrix.
			constexpr value_type determinant() const
			{
				return ((mA * mD) - (mC * mB));
			}

			// Computes the inverse of this matrix using the product returned by `determinant`.
			constexpr basic_matrix_2d inverse() const
			{
				return inverse(determinant());
			}
		
			// Computes the inverse of this matrix using the `determinant_value` specified.
			constexpr basic_matrix_2d inverse(value_type determinant_value) const
			{
				const auto a = (mD / determinant_value);
				const auto b = (-mB / determinant_value);
				const auto c = (-mC / determinant_value);
				const auto d = (mA / determinant_value);
		
				const auto tx = (((mC * mTy) - (mD * mTx)) / determinant_value);
				const auto ty = (((mB * mTx) - (mA * mTy)) / determinant_value);
		
				return { a, b, c, d, tx, ty };
			}
		
			// Resolves the X position of the point specified, using this matrix's geometric transformations.
			constexpr value_type transform_point_x(value_type point_x, value_type point_y) const
			{
				return ((mA * point_x) + (mC * point_y) + mTx);
			}

			// Resolves the Y position of the point specified, using this matrix's geometric transformations.
			constexpr value_type transform_point_y(value_type point_x, value_type point_y) const
			{
				return ((mB * point_x) + (mD * point_y) + mTy);
			}

			// Applies this matrix's geometric transformations to the specified point.
			constexpr std::pair<value_type, value_type> transform_point(value_type point_x, value_type point_y) const
			{
				return { transform_point_x(point_x, point_y), transform_point_y(point_x, point_y) };
			}

			// Applies this matrix's geometric transformations to the specified point.
			constexpr std::pair<value_type, value_type> transform_point(const std::pair<value_type, value_type>& point) const
			{
				return transform_point(point.first, point.second);
			}

			// Applies this matrix's geometric transformations to the specified point.
			constexpr typename std::enable_if<has_vector_type, VectorType>::type transform_point(const VectorType& point) const
			{
				auto result = transform_point(point.x, point.y);

				return VectorType { std::get<0>(result), std::get<1>(result) };
			}

			constexpr basic_matrix_2d multiply(const basic_matrix_2d& m) const
			{
				const auto a  = ((this->mA * m.mA) + (this->mB * m.mC));
				const auto b  = ((this->mA * m.mB) + (this->mB * m.mD));
				const auto c  = ((this->mC * m.mA) + (this->mD * m.mC));
				const auto d  = ((this->mC * m.mB) + (this->mD * m.mD));
				const auto tx = ((this->mTx * m.mA) + (this->mTy * m.mC));
				const auto ty = ((this->mTx * m.mB) + (this->mTy * m.mD));

				return basic_matrix_2d { a, b, c, d, tx, ty };
			}

			constexpr basic_matrix_2d transpose() const
			{
				return basic_matrix_2d
				{
					mA, mC, mTx,
					mB, mD, mTy
				};
			}

			// Accesses a matrix element by index.
#ifdef IMMUTABLEOCTET_MATRIX2D_TYPE_PUNNING
			constexpr value_type operator[](std::size_t index) const // const value_type&
			{
				assert(index < size());

				return raw_values[index];

                /*
                // Alternative implementation:
				switch (index)
				{
					case 0: return mA;
					case 1: return mB;
					case 2: return mC;
					case 3: return mD;
					case 4: return mTx;
					case 5: return mTy;
                    default: return value_type {};
				}
                */
			}
#endif

#ifdef IMMUTABLEOCTET_MATRIX2D_TYPE_PUNNING
			constexpr value_type& operator[](std::size_t index)
			{
				assert(index < size());

				return raw_values[index];

                /*
                // Alternative implementation:
				switch (index)
				{
					case 0: return mA;
					case 1: return mB;
					case 2: return mC;
					case 3: return mD;
					case 4: return mTx;
					case 5: return mTy;
				}
                */
			}
#endif

			// Creates a copy of this matrix then concatenates it with `m`, returning the result.
			constexpr basic_matrix_2d operator+(const basic_matrix_2d& m) const
			{
				return concatenate(m);
			}

			constexpr std::pair<value_type, value_type> operator*(const std::pair<value_type, value_type>& point) const
			{
				return transform_point(point);
			}

			constexpr typename std::enable_if<has_vector_type, VectorType>::type operator*(const VectorType& point) const
			{
				return transform_point(point);
			}

			constexpr basic_matrix_2d operator*(const basic_matrix_2d& m) const
			{
				return multiply(m);
			}
		protected:

#ifdef IMMUTABLEOCTET_MATRIX2D_TYPE_PUNNING
			union
			{
				struct
				{
					// Column  0    1    // Row
					value_type mA,  mB;  // 0
					value_type mC,  mD;  // 1
					value_type mTx, mTy; // 2
				};

				value_type raw_values[6]; // 3x2
			};
#else
			// Column  0    1    // Row
			value_type mA,  mB;  // 0
			value_type mC,  mD;  // 1
			value_type mTx, mTy; // 2
#endif
	};
}

#endif