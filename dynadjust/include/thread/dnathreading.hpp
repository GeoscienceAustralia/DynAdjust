//============================================================================
// Name         : dnathreading.hpp
// Author       : Roger Fraser
// Contributors :
// Version      : 1.00
// Copyright    : Copyright 2017 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//               
//                http ://www.apache.org/licenses/LICENSE-2.0
//               
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : DynAdjust multi-threading helps
//============================================================================

#ifndef DNATHREADING_H_
#define DNATHREADING_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD) 
#pragma message("  " __FILE__) 
#endif
#endif

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

using namespace std;
using namespace boost;

template<typename T>
class concurrent_block_adjustment
{
public:
	concurrent_block_adjustment() 
		: forward_is_running(false), reverse_is_running(false), combine_is_running(false) {
	}
	virtual ~concurrent_block_adjustment() {}

	inline void resize_forward_run(const T size) {
		// a lock is not needed since this function is called 
		// prior to creating threads 
		//lock_guard<mutex> lock_fwd(forward_mutex);
		v_fwd_blocks.resize(size);
		T id(0);
		for_each(v_fwd_blocks.begin(), v_fwd_blocks.end(),
			[this, &id] (sequential_adj_t<T, bool>& a) { 
				a.block_index = id++;
		});
	}

	inline void resize_reverse_run(const T size) {
		// a lock is not needed since this function is called 
		// prior to creating threads 
		//lock_guard<mutex> lock_rev(reverse_mutex);
		v_rev_blocks.resize(size);
		T id(0);
		for_each(v_rev_blocks.begin(), v_rev_blocks.end(),
			[this, &id] (sequential_adj_t<T, bool>& a) { 
				a.block_index = id++;
		});
	}
	inline void resize_runs(const T size) {
		resize_forward_run(size);
		resize_reverse_run(size);
	}

	inline void reset_forward_run() {
		// a lock is not needed since this function is called 
		// prior to creating threads 
		//lock_guard<mutex> lock_fwd(forward_mutex);
		for_each(v_fwd_blocks.begin(), v_fwd_blocks.end(), 
			boost::mem_fn(&sequential_adj::nosolution)); 
	}

	inline void reset_reverse_run() {
		// a lock is not needed since this function is called 
		// prior to creating threads 
		//lock_guard<mutex> lock_rev(reverse_mutex);
		for_each(v_rev_blocks.begin(), v_rev_blocks.end(), 
			boost::mem_fn(&sequential_adj::nosolution)); 
	}	

	inline void reset_adjustment_runs() {
		reset_forward_run();
		reset_reverse_run();
	}

	inline string print_adjusted_forward_blocks() {
		stringstream ss("");
		boost::lock_guard<boost::mutex> lock(forward_mutex);
		ss << "Forward blocks (" << v_fwd_blocks.size() << "):" << endl;
		for (UINT32 i=0; i<v_fwd_blocks.size(); ++i)
			if (v_fwd_blocks.at(i).adjusted)
				ss << v_fwd_blocks.at(i).block_index << " ";
		return ss.str();
	}

	inline string print_adjusted_reverse_blocks() {
		boost::lock_guard<boost::mutex> lock(reverse_mutex);
		stringstream ss;
		ss << "Reverse blocks (" << v_fwd_blocks.size() << "):" << endl;
		for (UINT32 i=0; i<v_rev_blocks.size(); ++i)
			if (v_rev_blocks.at(i).adjusted)
				ss << v_rev_blocks.at(i).block_index << " ";
		return ss.str();
	}

	inline string print_adjusted_blocks() {
		stringstream ss;
		ss << endl <<
			print_adjusted_forward_blocks() << endl <<
			print_adjusted_reverse_blocks() << endl;
		return ss.str();
	}

	// Set started state to true
	inline void combine_run_started() {					// combine run
		boost::lock_guard<boost::mutex> lock(cmb_state_mutex);
		combine_is_running = true; 
	}
	inline void forward_run_started() {					// forward run 
		boost::lock_guard<boost::mutex> lock(fwd_state_mutex);
		forward_is_running = true; 
	}
	inline void reverse_run_started() {					// reverse run
		boost::lock_guard<boost::mutex> lock(rev_state_mutex);
		reverse_is_running = true; 
	}
	
	// Set finished state to true
	inline void combine_run_stopped() {					// combine run
		boost::lock_guard<boost::mutex> lock(cmb_state_mutex);
		combine_is_running = false;						
	}													
	inline void forward_run_finished() {				// forward run
		boost::lock_guard<boost::mutex> lock(fwd_state_mutex);
		forward_is_running = false;						
	}													
	inline void reverse_run_finished() {				// reverse run
		boost::lock_guard<boost::mutex> lock(rev_state_mutex);
		reverse_is_running = false; 
	}
	
	// Get running state
	inline bool get_combine_state() {		// combine run
		boost::lock_guard<boost::mutex> lock(cmb_state_mutex);
		return combine_is_running;				  
	}											  
	inline bool get_forward_state() {		// forward run
		boost::lock_guard<boost::mutex> lock(fwd_state_mutex);
		return forward_is_running;				  
	}											  
	inline bool get_reverse_state() {		// reverse run
		boost::lock_guard<boost::mutex> lock(rev_state_mutex);
		return reverse_is_running; 
	}
	
	inline bool is_forward_finished() {		// forward run
		boost::lock_guard<boost::mutex> lock(fwd_state_mutex);
		return forward_is_running == false;				  
	}											  
	inline bool is_reverse_finished() {		// reverse run
		boost::lock_guard<boost::mutex> lock(rev_state_mutex);
		return reverse_is_running == false;				  
	}											  
	
	// It is dangerous to aquire a lock if one is already held 
	// (e.g. lock fwd_state_mutex, then lock rev_state_mutex)
	// so use unlock for forward before locking reverse
	inline bool any_still_running() { 
		if (get_forward_state())
			return true;		
		boost::lock_guard<boost::mutex> lock_r(rev_state_mutex);
		return reverse_is_running; 
	}
	inline bool finished_all_runs() { 
		if (get_forward_state())
			return false;		
		boost::lock_guard<boost::mutex> lock_r(rev_state_mutex);
		return !reverse_is_running; 
	}
	
	inline void set_forward_block_adjusted(const T& block) {
		boost::lock_guard<boost::mutex> lock(forward_mutex);
		v_fwd_blocks.at(block).solution();
	}
	inline void set_reverse_block_adjusted(const T& block) {
		boost::lock_guard<boost::mutex> lock(reverse_mutex);
		v_rev_blocks.at(block).solution();
	}
	
	inline bool forward_block_adjusted(const T& block) {
		boost::lock_guard<boost::mutex> lock(forward_mutex);
		return v_fwd_blocks.at(block).adjusted;
	}	
	inline bool reverse_block_adjusted(const T& block) {
		boost::lock_guard<boost::mutex> lock(reverse_mutex);
		return v_rev_blocks.at(block).adjusted;
	}
private:
	bool forward_is_running;
	bool reverse_is_running;
	bool combine_is_running;

	vector< sequential_adj_t<T, bool> > v_fwd_blocks;
	vector< sequential_adj_t<T, bool> > v_rev_blocks;

	mutable boost::mutex forward_mutex, reverse_mutex, combine_mutex;
	mutable boost::mutex fwd_state_mutex, rev_state_mutex, cmb_state_mutex;

};


template<typename T>
class concurrent_queue
{
public:
	concurrent_queue()
		: more_blocks_coming_(true) {}
	virtual ~concurrent_queue() {}

	inline void push_data(const vector<T>& data)
	{
		boost::lock_guard<boost::mutex> ql(queue_mutex);
		// push a vector of data elements
		for_each(
			data.begin(), 
			data.end(), 
			[this](const T& t){			// use lambda expression
				this->the_queue.push(t);
			}
		);
	}

	inline void push_and_notify(const T data)
	{
		boost::lock_guard<boost::mutex> ql(queue_mutex);
		the_queue.push(data);
		// notify all (potential) threads waiting on the_queue
		the_condition.notify_all();
	}

	inline bool not_empty() const
	{
		boost::lock_guard<boost::mutex> ql(queue_mutex);
		return !the_queue.empty();
	}

	inline bool is_empty() const
	{
		boost::lock_guard<boost::mutex> ql(queue_mutex);
		return the_queue.empty();
	}

	inline bool front_and_pop(T& popped_value)
	{
		boost::lock_guard<boost::mutex> ql(queue_mutex);
		if (the_queue.empty())
			return false;        
		popped_value = the_queue.front();
		the_queue.pop();
		return true;
	}

	inline void notify_all() 
	{
		the_condition.notify_all();
	}

	inline bool wait_front_and_pop(T& popped_value)
	{
		// wait if the queue is empty
		boost::unique_lock<boost::mutex> ql(queue_mutex);
		the_condition.wait(ql, [this] {return !this->the_queue.empty();});
		if (the_queue.empty())
			return false;
		popped_value = the_queue.front();
		the_queue.pop();
		return true;
	}

	inline void wait_if_queue_is_empty() 
	{
		boost::unique_lock<boost::mutex> ql(queue_mutex);
		the_condition.wait(ql, [this] { return !this->the_queue.empty();});
	}

	inline void queue_exhausted()
	{
		boost::lock_guard<boost::mutex> ql(exhausted_mutex);
		more_blocks_coming_ = false;
		the_condition.notify_all();
	}

	inline void reset_blocks_coming() 
	{ 
		boost::lock_guard<boost::mutex> el(exhausted_mutex);
		more_blocks_coming_ = true; 
	}

	inline bool more_blocks_coming() 
	{ 
		boost::lock_guard<boost::mutex> el(exhausted_mutex);
		return more_blocks_coming_; 
	}

	inline bool is_queue_exhausted()
	{
		// First, check that there are no more blocks coming
		if (more_blocks_coming())
			return false;

		// Second, is there any data in the queue?  If so, then 
		// that's the tell-tale sign we've got more work to do.
		return is_empty();
	}

	

private:
	mutable boost::mutex queue_mutex;
	mutable boost::mutex exhausted_mutex;
	std::queue<T> the_queue;
	boost::condition_variable the_condition;

	bool more_blocks_coming_;
};

template<typename T = string>
class concurrent_ofstream
{
public:
	concurrent_ofstream() {}
	virtual ~concurrent_ofstream() {}

	void wrtie(std::ofstream& ofs, const T& data)
	{
		boost::lock_guard<boost::mutex> lock(the_mutex);
		if(ofs.good())
			ofs << data;
	}
private:
	mutable boost::mutex the_mutex;
};

template<typename T>
class message_bank
{
public:
	message_bank() {}
	virtual ~message_bank() {}

	// Set started state to true
	inline bool messagebank_empty() {
		boost::lock_guard<boost::mutex> lock(message_mutex);
		return messages.empty(); 
	}
	
	inline T get_message(const UINT32& i) const {
		boost::lock_guard<boost::mutex> lock(message_mutex);
		if (i >= messages.size())
		{
			T msg(0);
			return msg;
		}
		return messages.at(i);
	}
	
	inline void add_message(const T& msg) {
		boost::lock_guard<boost::mutex> lock(message_mutex);
		messages.push_back(msg);
	}
	
	inline void clear_messages() {
		boost::lock_guard<boost::mutex> lock(message_mutex);
		messages.clear();
	}

private:
	vector<T> messages;
	mutable boost::mutex message_mutex;
};

template<typename T>
class protected_var
{
public:
	protected_var(const protected_var& v)
		:  protected_var(v, boost::mutex::scoped_lock(v.var_mutex)) {}

	protected_var() {}
	virtual ~protected_var() {}

	inline void write(const T& var) {
		boost::lock_guard<boost::mutex> lock(var_mutex);
		var_ = var;
	}
	
	inline T read() const {
		boost::lock_guard<boost::mutex> lock(var_mutex);
		return var_;
	}

private:
	protected_var(const protected_var& pv, const boost::mutex::scoped_lock& lock)
		:  var_(pv.var_) {}


private:
	T var_;
	mutable boost::mutex var_mutex;
};

#endif /* DNATHREADING_H_ */
