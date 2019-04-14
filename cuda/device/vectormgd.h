#pragma once

namespace ldpc
{


	template <typename T>
	class const_iterator_mgd
	{
	public:
		const_iterator_mgd(T *pContainer): mContainer(pContainer){}
		const_iterator_mgd operator++() { ++mContainer; return *this; }
		bool operator!=(const const_iterator_mgd& pOther) const { return mContainer != pOther.mContainer; }
		const T& operator*() const { return *mContainer; }
	private: 
		T* mContainer;
	};

	template <typename T>
	class iterator_mgd
	{
	public:
		iterator_mgd(T *pContainer): mContainer(pContainer){}
		iterator_mgd operator++() { ++mContainer; return *this; }
		bool operator!=(const iterator_mgd& pOther) const { return mContainer != pOther.mContainer; }
		T& operator*() const { return *mContainer; }
	private: 
		T* mContainer;
	};


	template <typename T>
	class vector_mgd
	{
	public:
		using iterator = iterator_mgd<T>;
		using const_iterator = const_iterator_mgd<T>;

		const_iterator begin() const { return const_iterator(mBuffer); }
		const_iterator end() const { return const_iterator(mBuffer + mLength); }
		iterator begin() { return iterator(mBuffer); }
		iterator end() { return iterator(mBuffer + mLength); }


		explicit vector_mgd(int pCap) : mCapacity(pCap), mLength(0)
		{
			if (cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity) != cudaSuccess)
			{
				throw; //TODO: throw exception
			}
			
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

		explicit vector_mgd(int pCap, const T& pVal) //init constructor
		: mCapacity(pCap)
		, mLength(0)
		{
			try
			{
				if (cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity) != cudaSuccess) { throw; }
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
				
				//continue exception
				throw;
			}
		}

		explicit vector_mgd(const vector_mgd& pCopy) //Copy constructor
		: mCapacity(pCopy.size())
		, mLength(0)
		{
			try
			{
				if (cudaMallocManaged(&mBuffer, sizeof(T)*mCapacity) != cudaSuccess) { throw; }
				for(size_t i = 0; i < pCopy.size(); ++i)
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
				
				//continue exception
				throw;
			}
		}

		inline vector_mgd& operator=(const vector_mgd& pCopy) { return *this; }
		inline const T& operator[](int pIndex) const {	return mBuffer[pIndex];	}
		inline T& operator[](int pIndex) { return mBuffer[pIndex]; }

		void push_back(const T& pVal) //no resize!!
		{
			resize_if_req();
			new(mBuffer + mLength) T(pVal); //copy buffer into new buffer
			++mLength;
		}
		void pop_back()
		{
			--mLength;
			mBuffer[mLength].~T(); //call destructor
		}

		T& at(int pIndex)
		{
			if (pIndex >= mLength) //check out of bound index
			{ throw; } else { return mBuffer[pIndex]; }
		}

		void resize(int pNewCap) 
		{
			//alloc new buffer with new size
			T* newBuff;
			if (cudaMallocManaged(&newBuff, sizeof(T)*pNewCap) != cudaSuccess) { throw; }
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

		__host__ __device__ inline size_t size() const { return mLength; }
	private:
		void resize_if_req() 
		{
			if (mLength == mCapacity) { resize(mCapacity+1); }
		}

		size_t mCapacity;
		size_t mLength;
		T* mBuffer;
	};

}
