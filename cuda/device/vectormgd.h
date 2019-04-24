#pragma once

#include <iostream>
#include <algorithm>

namespace ldpc
{
	template<typename T>
	void swap(T& one, T& two)
	{
		T tmp(one);
		one = two;
		two = tmp;
	}

	template <typename T>
	class cudamgd_ptr
	{
	public:
		__host__ __device__ cudamgd_ptr()
		: mContainer(nullptr),  mRefCount(nullptr) {}

		__host__ explicit cudamgd_ptr(T* pVal) //init constructor with pointer to obj, only on host
		: mContainer(pVal),  mRefCount(nullptr)
		{
			try
			{
				cudaError_t result = cudaMallocManaged(&mRefCount, sizeof(size_t));
				if (result != cudaSuccess)
				{
					throw std::runtime_error(cudaGetErrorString(result));
				}
				*mRefCount = 1;

				mem_prefetch();
			}
			catch(...)
			{
				throw;
			}
		}

		__host__ cudamgd_ptr(const T& pVal) //init constructor, only on host
		: mContainer(nullptr), mRefCount(nullptr)
		{
			try
			{
				cudaError_t result = cudaMallocManaged(&mContainer, sizeof(T));
				if (result != cudaSuccess)
				{
					throw std::runtime_error(cudaGetErrorString(result));
				}

				new(mContainer) T(pVal);

				result = cudaMallocManaged(&mRefCount, sizeof(size_t));
				if (result != cudaSuccess)
				{
					throw std::runtime_error(cudaGetErrorString(result));
				}
				*mRefCount = 1;

				mem_prefetch();
			}
			catch(...)
			{
				throw;
			}
		}

		__host__ __device__ cudamgd_ptr(const cudamgd_ptr& pCopy) //copy constructor
		: mContainer(pCopy.mContainer), mRefCount(pCopy.mRefCount) { (*mRefCount)++; }

		__host__ __device__ ~cudamgd_ptr()
		{
			if (mRefCount != nullptr)
			{
				if (--(*mRefCount) == 0) //only delete original pointer
				{
					mContainer->~T();
					cudaFree(mContainer);
					cudaFree(mRefCount);
				}
			}
			else if (mRefCount == nullptr)
			{
				if (mContainer != nullptr)
				{
					mContainer->~T();
					cudaFree(mContainer);
				}
			}
		}

		//copy/move assignment operator
		__host__ __device__ cudamgd_ptr& operator=(cudamgd_ptr pCopy) noexcept
		{
			swap(mContainer, pCopy.mContainer);
			swap(mRefCount, pCopy.mRefCount);

			return *this;
		}

		__host__ void mem_prefetch() //Prefetch, i.e. move the data to the gpu, to reduce latency
		{
			cudaDeviceSynchronize();

			int dev = -1;
			cudaGetDevice(&dev);
			cudaMemPrefetchAsync(mContainer, sizeof(T), dev, NULL);
			cudaMemPrefetchAsync(mRefCount, sizeof(size_t), dev, NULL);
		}

		__host__ __device__ T& operator*() { return *mContainer; }
		__host__ __device__ T* operator->() { return mContainer; }

		__host__ __device__ T* get() { return mContainer; }
	private:
		T* mContainer;
		size_t* mRefCount;
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

		__host__ vector_mgd() //default constructor
		: mCapacity(1), mLength(0), mBuffer(nullptr)
		{
			cudaError_t result = cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
			if (result != cudaSuccess)
			{
				throw std::runtime_error(cudaGetErrorString(result));
			}

		}

		__host__ vector_mgd(int pCap)
		: vector_mgd(pCap, T()) {}

		__host__ vector_mgd(int pCap, const T& pVal) //init constructor
		: mCapacity(pCap), mLength(0), mBuffer(nullptr)
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
				//mem_prefetch();
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
		: mCapacity(pCopy.mCapacity), mLength(0), mBuffer(nullptr)
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
				//prefetch when copy
				//mem_prefetch();
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
			if (mBuffer != nullptr) { cudaFree(mBuffer); }
		}

		//copy/move assignment operator
		//replaces both assignment operators
		__host__ vector_mgd& operator=(vector_mgd pCopy) noexcept //make a copy to eliminate const
		{
			swap(mCapacity, pCopy.mCapacity);
			swap(mLength, pCopy.mLength);
			swap(mBuffer, pCopy.mBuffer);

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

			//mem_prefetch();
		}

		//Prefetch, i.e. move the data to the gpu, to reduce latency
		__host__ void mem_prefetch()
		{
			cudaDeviceSynchronize();
			int dev = -1;
			cudaGetDevice(&dev);

			cudaError_t result = cudaMemPrefetchAsync(mBuffer, sizeof(T)*mCapacity, dev, NULL);
			if (result != cudaSuccess)
			{
				throw std::runtime_error(cudaGetErrorString(result));
			}
		}

		__host__ __device__ const_iterator begin() const { return const_iterator(mBuffer); }
		__host__ __device__ const_iterator end() const { return const_iterator(mBuffer + mLength); }
		__host__ __device__ iterator begin() { return iterator(mBuffer); }
		__host__ __device__ iterator end() { return iterator(mBuffer + mLength); }

		__host__ __device__ size_t size() const { return mLength; }
		__host__ __device__ T* data() const { return mBuffer; }

	private:
		__host__ __device__ void resize_if_req()
		{
			#ifdef __CUDA_ARCH__ //no resize on device
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
