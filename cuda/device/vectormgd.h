#pragma once

namespace ldpc
{


	template <typename T>
	class iterator
	{
	public:
		iterator(T *pContainer): mContainer(pContainer){}
		iterator operator++() { ++mContainer; return *this; }
		bool operator!=(const iterator& other) const { return mContainer != other.mContainer; }
		const T& operator*() const { return *mContainer; }
	private:
		T* mContainer;
	};



	template <typename T>
	class vector_mgd
	{
	public:

		iterator<T> begin() const { return iterator<T>(mBuffer); }
		iterator<T> end() const { return iterator<T>(mBuffer + mLength); }


		explicit vector_mgd(int pCap) : mCapacity(pCap), mLength(0)
		{
			cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
		}
		~vector_mgd()
		{
			//manually call destructor on elements in reverse order
			for(size_t i = 0; i < mLength; ++i)
			{
				mBuffer[mLength-1-i].~T();
			}
			cudaFree(mBuffer);
		}

		explicit vector_mgd(const vector_mgd& pCopy) //Copy constructor
		: mCapacity(pCopy.size())
		, mLength(0)
		{
			cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity);
			for(size_t i = 0; i < pCopy.size(); ++i)
			{
				push_back(pCopy[i]);
			}
		}

		inline vector_mgd& operator=(const vector_mgd& pCopy) { return *this; }
		const T& operator[](int index) const {
			return mBuffer[index];
		}
		T& operator[](int index) {
			return mBuffer[index];
		}

		void push_back(const T& pVal) //no resize!!
		{
			new(mBuffer + mLength) T(pVal); //copy buffer into new buffer
			++mLength;
		}
		void pop_back()
		{
			--mLength;
			mBuffer[mLength].~T(); // call destructor
		}
		T& at(int index);

		__host__ __device__ inline size_t size() const { return mLength; }
	private:
		size_t mCapacity;
		size_t mLength;
		T* mBuffer;
	};

}
