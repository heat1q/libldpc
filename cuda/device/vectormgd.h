#pragma once

namespace ldpc
{
	template <typename T>
	class cudamgd_ptr
	{
	public:
		__host__ cudamgd_ptr(const T& pVal) //init constructor, only on host
			: cudamgd_ptr(pVal, 1) {}

		__host__ cudamgd_ptr(const T& pVal, size_t pSize) //init constructor, only on host
			: mContainer(nullptr),  mIsRef(false)
		{
			size_t index = 0;
			try
			{
				cudaError_t result = cudaMallocManaged(&mContainer, sizeof(T)*pSize);
				if (result != cudaSuccess)
				{
					throw std::runtime_error(cudaGetErrorString(result));
				}
				for (; index < pSize; ++index) {
					new(mContainer + index) T(pVal);
				}
			}
			catch(...)
			{
				for(size_t i = 0; i < index; ++i) //destroy already constructed elements
				{
					mContainer[index-1-i].~T();
				}

				throw;
			}
		}

		__host__ __device__ cudamgd_ptr(const cudamgd_ptr& pCopy) //copy constructor
			: mContainer(pCopy.mContainer), mIsRef(true) {
			printf("Copy Pointer\n");}

		__host__ __device__ ~cudamgd_ptr()
		{
			#ifndef __CUDA_ARCH__ //only delete on host
				if (!mIsRef && mContainer != nullptr) //only delete original pointer
				{
					cudaFree(mContainer);
					printf("Delete Pointer\n");
				}
			#endif
		}

		__host__ __device__ cudamgd_ptr& operator=(const cudamgd_ptr& pCopy) //assignment operator
		{
			if (this != &pCopy) // Avoid self assignment
			{
				//if non referenced pointer, then delete current buffer
				if(!mIsRef)
				{
					cudaFree(mContainer);
				}

				//copy data & ref count
				mContainer = pCopy.mContainer;
				mIsRef = true;
			}
			return *this;
		}

		__host__ __device__ T& operator*() { return *mContainer; }
		__host__ __device__ T* operator->() { return mContainer; }

	private:
		T* mContainer;
		bool mIsRef;
	};


	template <typename T>
	class cuda_const_iterator
	{
	public:
		__host__ __device__ cuda_const_iterator(T *pContainer): mContainer(pContainer){}
		__host__ __device__ cuda_const_iterator operator++() { ++mContainer; return *this; }
		__host__ __device__ bool operator!=(const cuda_const_iterator& pOther) const { return mContainer != pOther.mContainer; }
		__host__ __device__ const T& operator*() const { return *mContainer; }
	private:
		T* mContainer;
	};

	template <typename T>
	class cuda_iterator
	{
	public:
		__host__ __device__ cuda_iterator(T *pContainer): mContainer(pContainer){}
		__host__ __device__ cuda_iterator operator++() { ++mContainer; return *this; }
		__host__ __device__ bool operator!=(const cuda_iterator& pOther) const { return mContainer != pOther.mContainer; }
		__host__ __device__ T& operator*() const { return *mContainer; }
	private:
		T* mContainer;
	};


	template <typename T>
	class vector_mgd
	{
	public:
		using iterator = cuda_iterator<T>;
		using const_iterator = cuda_const_iterator<T>;

		__host__ vector_mgd(int pCap)
			: mCapacity(pCap), mLength(0)
		{
			cudaError_t result = cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
			if (result != cudaSuccess)
			{
				throw std::runtime_error(cudaGetErrorString(result));
			}

		}

		__host__ vector_mgd(int pCap, const T& pVal) //init constructor
			: mCapacity(pCap), mLength(0)
		{
			try
			{
				cudaError_t result = cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
				if (result != cudaSuccess)
				{
					throw std::runtime_error(cudaGetErrorString(result));
				}
				for(size_t i = 0; i < mCapacity; ++i)
				{
					push_back(pVal);
				}
			}
			catch(...)
			{
				for(size_t i = 0; i < mLength; ++i) //destroy already constructed elements
				{
					mBuffer[mLength-1-i].~T();
				}

				//set length to zero to avoid cleaning elements in destructor
				mLength = 0;

				//continue exception
				throw;
			}
		}

		__host__ __device__ vector_mgd(const vector_mgd& pCopy) //copy constructor for device & host
			: mCapacity(pCopy.mCapacity), mLength(0)
		{
			#ifdef __CUDA_ARCH__
				//in device code, copy allocates memory for gpu
				cudaError_t result = cudaMalloc(&mBuffer, sizeof(T)*mCapacity);
				if (result != cudaSuccess)
				{
					printf("Error: vector_mgd: %s\n", cudaGetErrorString(result));
					exit(EXIT_FAILURE);
				}
				for(size_t i = 0; i < pCopy.mLength; ++i)
				{
					push_back(pCopy[i]);
				}
			#else
				try
				{
					cudaError_t result = cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
					if (result != cudaSuccess)
					{
						throw std::runtime_error(cudaGetErrorString(result));
					}
					for(size_t i = 0; i < pCopy.mLength; ++i)
					{
						push_back(pCopy[i]);
					}
				}
				catch(...)
				{
					for(size_t i = 0; i < mLength; ++i) //destroy already constructed elements
					{
						mBuffer[mLength-1-i].~T();
					}

					//set length to zero to avoid cleaning elements in destructor
					mLength = 0;

					//continue exception
					throw;
				}
			#endif
		}

		__host__ __device__ ~vector_mgd()
		{
			//manually call destructor on elements in reverse order
			for(size_t i = 0; i < mLength; ++i)
			{
				mBuffer[mLength-1-i].~T();
			}
			cudaFree(mBuffer);
		}

		//Operators
		__host__ __device__ vector_mgd& operator=(const vector_mgd& pCopy)
		{
			//TODO
			return *this;
		}
		__host__ __device__ const T& operator[](int pIndex) const {	return mBuffer[pIndex];	}
		__host__ __device__ T& operator[](int pIndex) { return mBuffer[pIndex]; }

		__host__ __device__ void push_back(const T& pVal)
		{
			resize_if_req();
			new(mBuffer + mLength) T(pVal); //copy buffer into new buffer
			++mLength;
		}

		__host__ void pop_back()
		{
			--mLength;
			mBuffer[mLength].~T(); //call destructor
		}

		__host__ T& at(int pIndex)
		{
			if (pIndex >= mLength) //check out of bound index
			{ throw std::runtime_error("bad index"); } else { return mBuffer[pIndex]; }
		}

		__host__ void resize(int pNewCap)
		{
			//alloc new buffer with new size
			T* newBuff;
			cudaError_t result = cudaMallocManaged(&newBuff, sizeof(T)*pNewCap);
			if (result != cudaSuccess)
			{
				throw std::runtime_error(cudaGetErrorString(result));
			}
			size_t newLen = 0;
			while (newLen < mLength && newLen < pNewCap) //check if new size is bigger or smaller than old
			{
				new(newBuff + newLen) T(mBuffer[newLen]);
				++newLen;
			}

			//destroy old elements
			for(size_t i = 0; i < mLength; ++i)
			{
				mBuffer[mLength-1-i].~T();
			}
			cudaFree(mBuffer);

			//assign new
			mBuffer = newBuff;
			mLength = newLen;
			mCapacity = pNewCap;
		}

		__host__ __device__ const_iterator begin() const { return const_iterator(mBuffer); }
		__host__ __device__ const_iterator end() const { return const_iterator(mBuffer + mLength); }
		__host__ __device__ iterator begin() { return iterator(mBuffer); }
		__host__ __device__ iterator end() { return iterator(mBuffer + mLength); }

		__host__ __device__ size_t size() const { return mLength; }

	private:
		__host__ __device__ void resize_if_req()
		{
			#ifdef __CUDA_ARCH__
				if (mLength == mCapacity)
				{
					printf("Error: vector_mgd: exceeds maximum capacity");
					exit(EXIT_FAILURE);
				}
			#else
				if (mLength == mCapacity) { resize(mCapacity+1); }
			#endif
		}

		size_t mCapacity;
		size_t mLength;
		T* mBuffer;
	};
}
