#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
int draw_unii(int i_min, int i_max, mt19937 &i_generator){
	//maybe just initialize once!
	//random_device i_rand_dev;
        //mt19937 i_generator(i_rand_dev());
	//until here
	uniform_int_distribution<int> i_distribution (i_min, i_max);
	return i_distribution(i_generator);
}

int make_weighted_draw(vector <double> &w_prob_vector, mt19937 &w_generator){
	//cout<<"Making weighted draw"<<endl;
	//maybe just initialize once!
	//random_device w_rand_dev;
        //mt19937 w_generator(w_rand_dev());
	//until here
	discrete_distribution<int> w_distribution (w_prob_vector.begin(), w_prob_vector.end());
	//cout<<"Done making weighted draw"<<endl;
	return w_distribution(w_generator);
}



class genotype_matrix{
	vector <char> genotype_array;
	vector <char> *ga_pointer;
	//only temp array is variable...
	vector <char> temp_array;
	vector <char> *ta_pointer;
	vector <int> position_array;
	vector <double> prob_vector;
	//position in bp
	int sd_site;
	//position in array
	int sd_position;
	int y_size;
	int x_size;
	int temp_x_size;
	int temp_last_position;
	int genome_size;
	//
	int m_rec_prob_window_size;
	//rec_prob vector:
	vector <double> mean_rec_prob_vector;
	//count the number of recombination vectors to get the mean of mean_rec_prob_vector
	int mean_rec_prob_vector_counter=0;
	//takes vector of sexes of individuals and adds them in this order to the genotype array...
	//determines left positions in for block var array and recombination probs
	vector <int> recomb_left_positions_var;		
	vector <double> recomb_probs_var;		
	public:
	genotype_matrix(vector <char> &ind_sexes, int gen_size, int sd_s, int mean_rec_prob_window_size){
		y_size=ind_sexes.size()*2;
		x_size=1;
		genome_size=gen_size;
		sd_site=sd_s;
		sd_position=0;
		position_array.push_back(sd_site);
		temp_last_position=0;
		temp_x_size=0;
		for(int i_ind=0; i_ind<ind_sexes.size()*2; i_ind++){
			if(ind_sexes[i_ind]=='m'){
				genotype_array.push_back('1');
				genotype_array.push_back('0');
			}
			else if(ind_sexes[i_ind]=='f'){
				genotype_array.push_back('0');
				genotype_array.push_back('0');
			}
		}
		genotype_array.reserve(10000);
		temp_array.reserve(10000);
		ga_pointer=&genotype_array;
		ta_pointer=&temp_array;
		prob_vector.reserve(gen_size-1);
		prob_vector.assign(gen_size-1,1);
		//number of windows in mean rec prob vector
		m_rec_prob_window_size=mean_rec_prob_window_size;
		int mean_rec_prob_vector_size=gen_size/mean_rec_prob_window_size;
		mean_rec_prob_vector.reserve(mean_rec_prob_vector_size);
		mean_rec_prob_vector.assign(mean_rec_prob_vector_size,0);

	}
	~genotype_matrix(){
		ga_pointer=nullptr;
		ta_pointer=nullptr;
	}
	//adds a single recombined gamete to the temp array, where ind_array.first defines the first and ind_array.second the second individual and left_rec_site the position (in position_array) of the first individual before the breakpoint. Array positions are 0based.
	void add_recombined_gamete(pair <int, int> ind_array_pos, const int &rg_window_size, const double &rg_min_ident, mt19937 &rc_randgen){
		//cout<<"Adding recombined gamete.."<<endl;
		make_rec_probs(ind_array_pos.first, ind_array_pos.second, rg_window_size, rg_min_ident);
		/*
		if(record_recprobs==true){
			//iterate over positions
			int curr_rec_wind=0;
			int curr_rec_wind_size=0;
			double curr_rec_wind_sum=0;
			for(int i_rp_v=0; i_rp_v<prob_vector.size(); i_rp_v++){
				if(curr_rec_wind_size<m_rec_prob_window_size){
					curr_rec_wind_sum+=prob_vector[i_rp_v];	
					curr_rec_wind_size++;
				}
				else{
					mean_rec_prob_vector[curr_rec_wind]+=curr_rec_wind_sum/curr_rec_wind_size;
					curr_rec_wind++;
					curr_rec_wind_size=0;
					curr_rec_wind_sum=prob_vector[i_rp_v];	
				}
			}
			mean_rec_prob_vector_counter++;
		}
		*/
		int left_rec_site=draw_rec_position( rc_randgen );
		//cout<<"Left_rec_site: "<<left_rec_site<<endl;
		//int left_rec_site=0;
		//check if break_point is before first or after last position in position array, in that case add random gamete from either parent...
		if(left_rec_site<position_array.front() || left_rec_site>=position_array.back()){
			if(draw_unii(0,1, rc_randgen)==0){
				non_recombined_gamete(ind_array_pos.first);
			}
			else{
				non_recombined_gamete(ind_array_pos.second);
			}
		}
		//find genotype position from position array...
		else{
			//cout<<"Doing binary search..."<<endl;
			int arr_pos_left=0;
			int arr_pos_right=position_array.size()-1;
			while(arr_pos_right-arr_pos_left>1){
				if(left_rec_site==position_array[arr_pos_left]){
					break;
				}
				else if(left_rec_site==position_array[arr_pos_right]){
					arr_pos_left=arr_pos_right;
					break;
				}
				else if(left_rec_site>=position_array[arr_pos_left+((arr_pos_right-arr_pos_left)/2)]){
					arr_pos_left+=(arr_pos_right-arr_pos_left)/2;
				}
				else if(left_rec_site<position_array[arr_pos_left+((arr_pos_right-arr_pos_left)/2)]){
					arr_pos_right-=(arr_pos_right-arr_pos_left)/2;
				}
			//	cout<<"Rec position:"<<left_rec_site<<endl;
		//		cout<<"Rec_pos: "<<left_rec_site<<" arr_pos_left: "<<arr_pos_left<<", value: "<<position_array[arr_pos_left]<<" arr_pos_right: "<<arr_pos_right<<", value:"<<position_array[arr_pos_right]<<endl;
		//		print_position_array();
			}

			recombined_gamete(arr_pos_left ,ind_array_pos);
		}
		//cout<<"Done adding recombined gamete.."<<endl;
	}

	void add_recombined_gamete_random(pair <int, int> ind_array_pos_rand, mt19937 &rc_rand_randgen){
		int rand_lr_site=draw_unii(0, position_array.size()-1, rc_rand_randgen);
		if(rand_lr_site<position_array.front() || rand_lr_site>=position_array.back()){
			if(draw_unii(0,1, rc_rand_randgen)==0){
				non_recombined_gamete(ind_array_pos_rand.first);
			}
			else{
				non_recombined_gamete(ind_array_pos_rand.second);
			}
		}
		else{
			recombined_gamete(rand_lr_site, ind_array_pos_rand);
		}
	}


	void add_recombined_gamete_variable_sites(pair <int, int> ind_array_pos_var, const int &rg_window_size_var, const double &rg_min_ident_var, mt19937 &rc_randgen_var, double general_rec_prob){


		//TODO: do a random draw to test if recombination will take place
		vector <double> gr_probs={general_rec_prob, 1-general_rec_prob};
		if(make_weighted_draw(gr_probs, rc_randgen_var)==0){

			if(recomb_left_positions_var.size()==0){
				recomb_left_positions_var.emplace_back(-1);	
			}
			else{
				recomb_left_positions_var[0]=-1;
			}
			//cout<<"Adding gamete var sites..."<<endl;
			//1. determine the number of possible recombination points
			int rc_points=0;
			//where does the current block begin
			int curr_block_begin=0;	
			for(int i_rc_pos=0; i_rc_pos<position_array.size(); i_rc_pos++){
				if((*ga_pointer)[(ind_array_pos_var.first*x_size)+i_rc_pos] != (*ga_pointer)[(ind_array_pos_var.second*x_size)+i_rc_pos]){
	//				cout<<"Position "<<i_rc_pos<<" in Sequence "<<ind_array_pos_var.first<<" and Sequence "<<ind_array_pos_var.second<<" are variable."<<endl;
					double raw_prob=double(position_array[i_rc_pos])-curr_block_begin;
					if(recomb_probs_var.size()<rc_points+1){
						recomb_probs_var.emplace_back(raw_prob);
					}	
					else{
						recomb_probs_var[rc_points]=raw_prob;
					}
					curr_block_begin=position_array[i_rc_pos];
					rc_points++;
					if(recomb_left_positions_var.size()>rc_points){
						recomb_left_positions_var[rc_points]=i_rc_pos;
					}
					else{
						recomb_left_positions_var.emplace_back(i_rc_pos);
					}
				}
			}
			//add probability for last window as well...
			if(rc_points>0){
				double last_raw_prob=double(genome_size)-curr_block_begin-1;
				if(recomb_probs_var.size()<=rc_points){
					recomb_probs_var.emplace_back(last_raw_prob);
				}	
				else{
					recomb_probs_var[rc_points]=last_raw_prob;
				}
				//if necessary shrink both vectors to correct size
				while(recomb_probs_var.size()>rc_points+1){
					recomb_left_positions_var.pop_back();
					recomb_probs_var.pop_back();
				}
				//print rec prob vector
				/*
				cout<<"Made rec prob vector."<<endl;
				double rec_prob_sum=0;
				for(int i_v_pp=0; i_v_pp<recomb_probs_var.size(); i_v_pp++){
					cout<<recomb_probs_var[i_v_pp]<<";";
					rec_prob_sum+=recomb_probs_var[i_v_pp];
				}

				cout<<endl;
				cout<<"Recprob sum: "<<rec_prob_sum<<",N recprobs:"<<recomb_probs_var.size()<<endl;
				for(int i_v_ps=0; i_v_ps<recomb_left_positions_var.size(); i_v_ps++){
					cout<<recomb_left_positions_var[i_v_ps]<<";";
				}
				cout<<endl;
				cout<<"N recsites: "<<recomb_left_positions_var.size()<<endl;

				*/
				//
				//TODO: iterate over recomb_probs_var and reduce rec probs for each block			
				for(int i_rc_red=0; i_rc_red<recomb_probs_var.size(); i_rc_red++){
					int block_left;
					if(recomb_left_positions_var[i_rc_red]!=-1){
						block_left=position_array[recomb_left_positions_var[i_rc_red]];
					}
					else{
						block_left=0;
					}
					int block_right=block_left+int(recomb_probs_var[i_rc_red]);
					//iterate to left of block
					if(block_left>0){
						int l_red_it=i_rc_red;
						while(recomb_left_positions_var[l_red_it]!=-1 && position_array[recomb_left_positions_var[l_red_it]]>block_left-rg_window_size_var){
							//cout<<"Found left var site..."<<endl;
							double l_red;	
							int l_wind_begin=position_array[recomb_left_positions_var[l_red_it]];
							//cout<<"L wind begin:"<<l_wind_begin<<endl;
							int l_wind_end=l_wind_begin+rg_window_size_var;
							//cout<<"L wind end:"<<l_wind_end<<endl;
							//calculate reduction based on window size and overlap
							//if block is smaller than window end reduce proportionate to block size
							if(l_wind_end>block_right){
								l_red=recomb_probs_var[i_rc_red]*1/double(rg_window_size_var);	
							}
							//otherwise reduce proportionate to overlap with block
							else{
								l_red=(l_wind_end-block_left)*1/double(rg_window_size_var);
							}
							//cout<<"l_red:"<<l_red<<endl;
							recomb_probs_var[i_rc_red]-=l_red;
							l_red_it--;
						}
					}
					//iterate to right of block
					if(block_right<genome_size){
						int r_red_it=i_rc_red+1;
						while(r_red_it<recomb_left_positions_var.size() && position_array[recomb_left_positions_var[r_red_it]]<block_right+rg_window_size_var){
							//cout<<"Found right var site..."<<endl;
							double r_red;
							int r_wind_end=position_array[recomb_left_positions_var[r_red_it]];
							int r_wind_begin=r_wind_end-rg_window_size_var;
							if(r_wind_begin<block_left){
								r_red=recomb_probs_var[i_rc_red]*1/double(rg_window_size_var);
							}
							else{
								r_red=(block_right-r_wind_begin)*1/double(rg_window_size_var);
							}
							//cout<<"r_red:"<<r_red<<endl;
							recomb_probs_var[i_rc_red]-=r_red;
							r_red_it++;
						}
					}
					//TODO: check if window probs are below minimum threshold
					double window_min_threshold;
					//if the block is larger than two times the window it cannot completely be reduced to zero...
					if(block_right-block_left>2*rg_window_size_var){
						window_min_threshold=2*rg_window_size_var*double(rg_min_ident_var);
						if(recomb_probs_var[i_rc_red]<window_min_threshold){
							recomb_probs_var[i_rc_red]=block_right-block_left-(2*rg_window_size_var);
						}	
					}
					else{
						window_min_threshold=(block_right-block_left)*double(rg_min_ident_var);
						if(recomb_probs_var[i_rc_red]<window_min_threshold){
							recomb_probs_var[i_rc_red]=0;
						}	
						
					}


				}
				int rc_p_index=make_weighted_draw(recomb_probs_var, rc_randgen_var);
				if(rc_p_index==recomb_left_positions_var.front() || rc_p_index==recomb_left_positions_var.back()){
					if(draw_unii(0,1, rc_randgen_var)==0){
						non_recombined_gamete(ind_array_pos_var.first);
					}
					else{
						non_recombined_gamete(ind_array_pos_var.second);
					}
				}
				//3. make recombined gamete
				else{
					recombined_gamete(rc_p_index, ind_array_pos_var);
				}
			
			}
		/*
				//print rec prob vector
				cout<<"Change rec probs:"<<endl;
				rec_prob_sum=0;
				for(int i_v_pp=0; i_v_pp<recomb_probs_var.size(); i_v_pp++){
					cout<<recomb_probs_var[i_v_pp]<<";";
					rec_prob_sum+=recomb_probs_var[i_v_pp];
				}
				cout<<endl;
				*/


			//if both sequences are identical
			else{
					non_recombined_gamete(ind_array_pos_var.first);
			}

		}
		else{
				if(draw_unii(0,1, rc_randgen_var)==0){
					non_recombined_gamete(ind_array_pos_var.first);
				}
				else{
					non_recombined_gamete(ind_array_pos_var.second);
				}
		}
	}


	void print_position_array(){
		cout<<"Position array ("<<position_array.size()<<" bp):"<<endl;
		for(int i_p_arr=0; i_p_arr<position_array.size(); i_p_arr++){
				cout<<position_array[i_p_arr]<<";";
		}
		cout<<endl;
	}
	//adds mutations to the temp array, operates on position array
	void temp_add_mutations(double mutation_rate, double male_bias, vector <int> *mut_females, vector <int> *mut_males, mt19937 &am_randgen){
		//draw total number of mutations in males (use temp female array)
		int male_muts=mut_males->size()*2*genome_size*mutation_rate*male_bias;
		//draw total number of mutations in females
		int female_muts=mut_females->size()*2*genome_size*mutation_rate;
		//draw number of mutations in males and in females (use temp male array)
		int f_mut_count=0;
		int m_mut_count=0;
		while(f_mut_count+m_mut_count<male_muts+female_muts){
			//draw mutation position
			bool mutation_added=false;
			while(mutation_added==false){
				int mut_pos=draw_unii(0,genome_size-1, am_randgen);
				if(mut_pos>position_array.back()){
					int mut_ind_end;
					//cout<<"Adding mutation at end of array..."<<endl;
					if(f_mut_count<female_muts){
						//draw individual
						mut_ind_end=(*mut_females)[draw_unii(0, mut_females->size()-1, am_randgen)];
						f_mut_count++;
						}
					else{
						mut_ind_end=(*mut_males)[draw_unii(0, mut_males->size()-1, am_randgen)];
						m_mut_count++;
					}	
					//draw haplotype
					int mut_seq_end=mut_ind_end+draw_unii(0,1, am_randgen);
					position_array.push_back(mut_pos);
					add_site(temp_x_size, mut_seq_end);
					mutation_added=true;
				}
				else if(mut_pos<position_array.front()){
					int mut_ind_beg;
					//cout<<"Adding mutation at beginning of array..."<<endl;
					if(f_mut_count<female_muts){
						mut_ind_beg=(*mut_females)[draw_unii(0, mut_females->size()-1, am_randgen)];
						f_mut_count++;
					}
					else{
						mut_ind_beg=(*mut_males)[draw_unii(0, mut_males->size()-1, am_randgen)];
						m_mut_count++;
					}
					//int mut_ind_beg=draw_unii(0, y_size-1, am_randgen);
					position_array.insert(position_array.begin(), mut_pos);
					int mut_seq_beg=mut_ind_beg+draw_unii(0,1, am_randgen);
					add_site(0, mut_seq_beg);
					mutation_added=true;
					sd_position++;
				}
				else{
					int left_index=0;
					int right_index=position_array.size()-1;
					bool position_in_array=false;
					while(right_index-left_index>1){
						if((mut_pos==position_array[left_index]) || (mut_pos==position_array[right_index])){
							position_in_array=true;
							break;
						}
						int mid_index;
						if(right_index-left_index%2==0){
							mid_index=left_index+((right_index-left_index)/2);
						}
						else{
							mid_index=left_index+((right_index-left_index+1)/2);
						}
						/*
						if(right_index-left_index==2){
							left_index++;
						}
						*/
						if(mut_pos==position_array[mid_index]){
							position_in_array=true;
							break;
						}
						else if(mut_pos>position_array[mid_index]){
							left_index=mid_index;
						}
						else if(mut_pos<position_array[mid_index]){
							right_index=mid_index;
						}
						//cout<<"Mutation position: "<<mut_pos<<endl;
						//cout<<"Left index: "<<left_index<<", Mut pos: "<<position_array[left_index]<<endl;
						//cout<<"Right index: "<<right_index<<", Mut pos: "<<position_array[right_index]<<endl;
					}
					if(position_in_array==false){
						
						int mut_ind_mid;
						if(f_mut_count<female_muts){
						//draw individual
							mut_ind_mid=(*mut_females)[draw_unii(0, mut_females->size()-1, am_randgen)];
							f_mut_count++;
						}
						else{
							mut_ind_mid=(*mut_males)[draw_unii(0, mut_males->size()-1, am_randgen)];
							m_mut_count++;
						}	
						int mut_ind_seq=mut_ind_mid+draw_unii(0,1, am_randgen);
						position_array.insert(position_array.begin()+right_index, mut_pos);
						//cout<<"Added site in position array.."<<endl;
						add_site(right_index, mut_ind_seq);
						//cout<<"Added site in genotype matrix"<<endl;
						mutation_added=true;
						if(mut_pos<sd_site){
							sd_position++;
						}
					}
				}
			}
		}
	}
	void print_temp_array(){
		cout<<"Temp array(temp_x_size: "<<temp_x_size<<", vector size: "<<ta_pointer->size()<<", capacity: "<<ta_pointer->capacity()<<", temp last position: "<<temp_last_position<<")"<<endl;

		for(int i_t_arr=0; i_t_arr<temp_last_position; i_t_arr++){
			if((i_t_arr+1)%(temp_x_size)==0){
				cout<<(*ta_pointer)[i_t_arr]<<" "<<(i_t_arr)/temp_x_size<<endl;
			}
			else if(i_t_arr%sd_position==0 && (*ta_pointer)[i_t_arr]=='1'){
				cout<<"\033[0;31m"<<(*ta_pointer)[i_t_arr]<<"\033[0m";
			
			}
			else{
				cout<<(*ta_pointer)[i_t_arr];
			}
		}
		cout<<endl;
		cout<<"SD position: "<<sd_position<<endl;
	}
	//removes fixed sites from the temp array
	void temp_remove_fixed_sites(){
		//cout<<"Removing fixed sites.."<<endl;
		//iterate over sites in genotype array
		int i_f_site=0;
		while(i_f_site<temp_x_size){
		//	cout<<"Site: "<<i_f_site<<endl;
			//iterate over sequences
			char curr_nuc='n';
			bool fixed=true;
			for(int i_f_seq=i_f_site; i_f_seq<y_size*temp_x_size; i_f_seq+=temp_x_size){
				if(curr_nuc=='n'){
					curr_nuc=(*ta_pointer)[i_f_seq];
				}
				else if((*ta_pointer)[i_f_seq]!=curr_nuc){
					fixed=false;
					break;
				}
			}
			//remove fixed site and change iterator accordingly...
			if(fixed==true){
			//	cout<<"Removing fixed site at position "<<i_f_site<<endl;
				if(position_array[i_f_site]<sd_site){
					sd_position--;
				}
				remove_site(i_f_site);
				position_array.erase(position_array.begin()+i_f_site);
				//print_temp_array();
			}
			else{
				i_f_site++;
			}
		}
		//cout<<"Done removing fixed sites..."<<endl;
	}
	pair <char, char> return_genotype(int g_pos, pair <int, int> g_array_pos){
		return pair <char, char> ((*ga_pointer)[g_array_pos.first], (*ga_pointer)[g_array_pos.second] );
	}
	void swap_arrays(){
		x_size=temp_x_size;
		temp_x_size=0;
		temp_last_position=0;
		//TODO: play around with vector capacities...
		//genotype_array.swap(temp_array);
		swap(ga_pointer, ta_pointer);

	}
	void print_array(){
		cout<<"Genotype array:"<<", vector size: "<<ga_pointer->size()<<", capacity: "<<ga_pointer->capacity()<<")"<<endl;
		for(int i_arr=0; i_arr<y_size*x_size; i_arr++){
			if((i_arr+1)%(x_size)==0){
				cout<<(*ga_pointer)[i_arr]<<" "<<(i_arr)/x_size<<endl;
			}
			/*else if((i_arr)%(x_size+sd_position)==0){
				cout<<"\033[1;31m"<<(*ga_pointer)[i_arr]<<"\033[0m";
			}
			*/
			else{
				cout<<(*ga_pointer)[i_arr];
			}
		}
		cout<<endl;
	}
	void print_temp_vector(){
		cout<<"Temp x: "<<temp_x_size<<endl;
		cout<<"Temp vector ("<<ta_pointer->size()<<" sites)"<<endl;
		for(int i_t_vec=0; i_t_vec<ta_pointer->size(); i_t_vec++){
			cout<<(*ta_pointer)[i_t_vec]<<";";
		}
		cout<<endl;

	}
	char determine_sex(pair <int, int> sd_sequence){
		if((*ta_pointer)[(sd_sequence.first*temp_x_size)+sd_position]=='1' || (*ta_pointer)[(sd_sequence.second*temp_x_size)+sd_position]=='1'){
			return 'm';
		}
		else{
			return 'f';
		}
	}
	void make_rec_probs(int seq1, int seq2, const int &window_size, const double &min_ident){
		prob_vector.assign(prob_vector.size(),1);
		//iterate over position array
		for(int i_pa=0; i_pa<position_array.size(); i_pa++){
			//check if both sequences are different!! in position
			if(genotype_array[(x_size*seq1)+i_pa]!=genotype_array[(x_size*seq2)+i_pa]){
				//change prob vector here
				int decrease_begin;
				int decrease_end;
				int decrease_size=window_size;
				if(position_array[i_pa]-(window_size/2)<0){
					decrease_begin=0;
					decrease_size+=position_array[i_pa]-(window_size/2);					
				}
				else if(position_array[i_pa]+(window_size/2)>prob_vector.size()){
					decrease_end=prob_vector.size();
					decrease_size-=(position_array[i_pa]+window_size)-prob_vector.size();					
				}
				else{
					decrease_begin=position_array[i_pa]-(decrease_size/2);
					decrease_end=position_array[i_pa]+(decrease_size/2);
				}
				//iterate over window
				for(int i_dec=decrease_begin; i_dec<=decrease_end; i_dec++){
					//this model: rec probability relative to seq identity in window
					prob_vector[i_dec]-=1/double(decrease_size);
					if(prob_vector[i_dec]<min_ident){
						prob_vector[i_dec]=0;
					}
				}
			}		
		}
	}
	int draw_rec_position(mt19937 &rd_rand_gen){
			return make_weighted_draw(prob_vector, rd_rand_gen);
	}
	void write_size_of_sdr(vector <int> *sdr_males, ofstream &sdr_file){
		//move left from position array and determine if mutations are fixed in males
		bool leftmost_fixed=true;
		int left_array_position=sd_position;
		while(leftmost_fixed==true){
			left_array_position--;
			if (left_array_position>=0){
				//iterate over males
				for(int ls_mi=0; ls_mi<sdr_males->size(); ls_mi++){
					char lm_gen1=(*ga_pointer)[((*sdr_males)[ls_mi]*x_size)+left_array_position];
					char lm_gen2=(*ga_pointer)[(((*sdr_males)[ls_mi]+1)*x_size)+left_array_position];
					if (lm_gen1==lm_gen2){
						leftmost_fixed=false;
						break;
					}
				}
			}
			else{
				leftmost_fixed=false;
			}
		}
		left_array_position++;
		bool rightmost_fixed=true;
		int right_array_position=sd_position;
		while(rightmost_fixed==true){
			right_array_position++;
			//get size of genotype array
			if (right_array_position<=x_size-1){
				//iterate over males
				for(int rs_mi=0; rs_mi<sdr_males->size(); rs_mi++){
					char rm_gen1=(*ga_pointer)[((*sdr_males)[rs_mi]*x_size)+right_array_position];
					char rm_gen2=(*ga_pointer)[(((*sdr_males)[rs_mi]+1)*x_size)+right_array_position];
					if (rm_gen1==rm_gen2){
						rightmost_fixed=false;
						break;
					}
				}
			}
			else{
				rightmost_fixed=false;
			}
		}
		left_array_position++;
		int sdr_size=position_array[right_array_position]-position_array[left_array_position];
		sdr_file<<"SDR size:"<<sdr_size<<endl;
	}
	
	void write_rec_probs(){
		//estimate mean for each window
		vector <double> rp_mean(mean_rec_prob_vector.size(), 0);				
		//write estimate to output file
		for(int i_rp_m=0; i_rp_m<mean_rec_prob_vector.size(); i_rp_m++){
			rp_mean[i_rp_m]=mean_rec_prob_vector[i_rp_m]/mean_rec_prob_vector_counter;
		}
		//empty vector again afterwards
		mean_rec_prob_vector.assign(mean_rec_prob_vector.size(),0);
		mean_rec_prob_vector_counter=0;
	}	
	//m/f FST
	void write_window_fst(int fst_window_size, vector <int> *fst_m_pointer, vector <int> *fst_f_pointer, const int &fst_gen, ostream &fst_outfile){
		//make vector of windows
		int n_windows=genome_size/fst_window_size;
		//vector <double> fst_per_window(n_windows, 0);
		//just cout for now, later into an output file...
	//	cout<<"FST: ";
		//iterate over position array
		fst_outfile<<fst_gen<<"\t"<<position_array.size();
		

		int curr_window_start=0;
		int curr_window_end=fst_window_size;
		int fst_window_counter=0;
		int curr_window_sites=0;
		double curr_window_fst=0;
		for(int i_fst=0; i_fst<position_array.size(); i_fst++){
			int n_m=fst_m_pointer->size();
			int n_f=fst_f_pointer->size();
			//iterate over males
			int m_all_count=0;			
			for(int i_m_a=0; i_m_a<n_m; i_m_a++){
				if((*ga_pointer)[((*fst_m_pointer)[i_m_a]*x_size)+i_fst]=='1'){
					m_all_count++;
				}
				if((*ga_pointer)[(((*fst_m_pointer)[i_m_a]+1)*x_size)+i_fst]=='1'){
					m_all_count++;
				}
			}
			double m_a_freq=double(m_all_count)/(2*n_m);		

			double Hm=1-((m_a_freq*m_a_freq)+((1-m_a_freq)*(1-m_a_freq)));
			//iterate over females
			int f_all_count=0;
			for(int i_f_a=0; i_f_a<n_f; i_f_a++){
				if((*ga_pointer)[((*fst_f_pointer)[i_f_a]*x_size)+i_fst]=='1'){
					f_all_count++;
				}
				if((*ga_pointer)[(((*fst_f_pointer)[i_f_a]+1)*x_size)+i_fst]=='1'){
					f_all_count++;
				}
			}
			double f_a_freq=double(f_all_count)/(2*n_f);		
			double Hf=1-((f_a_freq*f_a_freq)+((1-f_a_freq)*(1-f_a_freq)));
			//estimate expected heterozygosity over all samples
			double comb_a_freq=double(f_all_count+m_all_count)/(2*(n_f+n_m));
			double HT=1-((comb_a_freq*comb_a_freq)+((1-comb_a_freq)*(1-comb_a_freq)));
			double site_fst=(HT-((n_m*Hm)+(n_f*Hf))/(n_m+n_f))/HT;
			if(position_array[i_fst]<curr_window_end){
				curr_window_fst+=site_fst;
				curr_window_sites++;
			}
			else{	
				curr_window_start=curr_window_end;
				curr_window_end+=fst_window_size;
				if(curr_window_sites>0){
				//	fst_per_window[fst_window_counter]=curr_window_fst/curr_window_sites;
					fst_outfile<<"\t"<<curr_window_fst/double(curr_window_sites);
				}
				else{
				//	fst_per_window[fst_window_counter]=0;
					fst_outfile<<"\tNA";
				}
				curr_window_fst=0;
				curr_window_sites=0;
				fst_window_counter++;
				i_fst--;
			}
		}
		//write last window
		if(curr_window_sites>0){
			//fst_per_window[fst_window_counter]=curr_window_fst/curr_window_sites;
			fst_outfile<<"\t"<<curr_window_fst/double(curr_window_sites)<<endl;
		}
		else{
			//fst_per_window[fst_window_counter]=0;
			fst_outfile<<"\t"<<"NA"<<endl;
		}
	}
	void write_window_dxy(int dxy_window_size, vector <int> *dxx_f_pointer, vector <int> *dxy_m_pointer, ostream &dxx_outfile, ostream &dxy_outfile, int dxy_gen){
		int curr_window_end=dxy_window_size;
		int n_windows=genome_size/dxy_window_size;
		dxy_outfile<<dxy_gen;
		dxx_outfile<<dxy_gen;
		int dxy_pos_array_count=0;
		//iterate over windows
		for(int i_dxy_window=0; i_dxy_window<n_windows; i_dxy_window++){
			double curr_window_dxy=0;
			double curr_window_dxx=0;
			while(position_array[dxy_pos_array_count]<curr_window_end){
		
				//iterate over males
				for(int i_dxy_m=0; i_dxy_m<dxy_m_pointer->size(); i_dxy_m++){
					if( (*ga_pointer)[ ((*dxy_m_pointer)[i_dxy_m]*x_size)+dxy_pos_array_count]!= (*ga_pointer)[(((*dxy_m_pointer)[i_dxy_m]+1)*x_size)+dxy_pos_array_count] ){
					curr_window_dxy+=1/(double(dxy_window_size)*dxy_m_pointer->size());
					}
				}
				//iterate over females
				for(int i_dxx_f=0; i_dxx_f<dxx_f_pointer->size(); i_dxx_f++){
					if( (*ga_pointer)[ ((*dxx_f_pointer)[i_dxx_f]*x_size)+dxy_pos_array_count]!= (*ga_pointer)[(((*dxx_f_pointer)[i_dxx_f]+1)*x_size)+dxy_pos_array_count] ){
					curr_window_dxx+=1/(double(dxy_window_size)*dxx_f_pointer->size());
					}
				}
				dxy_pos_array_count++;
			}
			dxy_outfile<<"\t"<<curr_window_dxy;
			dxx_outfile<<"\t"<<curr_window_dxx;
			curr_window_end+=dxy_window_size;
		}
		dxy_outfile<<endl;
		dxx_outfile<<endl;
	}

	private:
		//fill the temp array by overwriting existing genotypes and only change the size when we go past the end and in the second step, when mutations are added and fixed sites are removed!
		void recombined_gamete(int ga_rec_pos_left, pair <int, int> &seq1seq2ids){
		//	cout<<"Adding rec gamete"<<endl;
			//first find out if the temp array is large enough, if not push back zeros until it is.
			//cout<<"Adding rec gamete! ta_vector size: "<<ta_pointer->size()<<", temp_last_position: "<<temp_last_position<<endl;
			while(ta_pointer->size()<temp_last_position+x_size){
				ta_pointer->emplace_back('0');
			}
			//positions in genotype array
			int seq1_start_position=(seq1seq2ids.first*x_size);
			int seq2_start_position=(seq1seq2ids.second*x_size)+ga_rec_pos_left+1;

			for(int i_rec_pos=seq1_start_position; i_rec_pos<=seq1_start_position+ga_rec_pos_left; i_rec_pos++){
				(*ta_pointer)[temp_last_position]=(*ga_pointer)[i_rec_pos];
				temp_last_position++;
			}
			//for(int i_rec_pos_2=ga_rec_pos_left+1; i_rec_pos_2<x_size; i_rec_pos_2++){
			for(int i_rec_pos_2=seq2_start_position; i_rec_pos_2<(seq1seq2ids.second+1)*x_size; i_rec_pos_2++){
				(*ta_pointer)[temp_last_position]=(*ga_pointer)[i_rec_pos_2];
				temp_last_position++;
			}
			temp_x_size=x_size;
		//	cout<<"Done adding rec gamete"<<endl;
		}
		void non_recombined_gamete(int seq_id){
			//cout<<"Adding nonrec gamete! ta_vector size: "<<ta_pointer->size()<<", temp_last_position: "<<temp_last_position<<endl;
			int seq_start_position=seq_id*x_size;
			while(ta_pointer->size()<temp_last_position+x_size){
				ta_pointer->emplace_back('0');
			//	cout<<"Emplacing zero to temp vector, new temp vector size: "<<ta_pointer->size()<<endl;
			}


			for(int i_non_rec_pos=seq_start_position; i_non_rec_pos<(seq_id+1)*x_size; i_non_rec_pos++){
			//	cout<<"Copying site "<<i_non_rec_pos<<" from genotype array position "<<i_non_rec_pos<<" to temp array position "<<temp_last_position<<endl;
				(*ta_pointer)[temp_last_position]=(*ga_pointer)[i_non_rec_pos];
				temp_last_position++;
			//	cout<<"Temp last position: "<<temp_last_position<<endl;
			}
			temp_x_size=x_size;
		}
		//here position is position in temp array!
		void add_site(int position, int seq_ID){
		//	cout<<"Adding site at position "<<position<<" in sequence "<<seq_ID<<" out of a total of "<<temp_x_size<<" positions in a vector of size "<<ta_pointer->size()<<" which ends at position "<<temp_last_position<<endl;
			//iterate over n_sequences
			int curr_pos=position;
			for(int i_add_seq=0; i_add_seq<y_size; i_add_seq++){
				if(seq_ID==i_add_seq){
					ta_pointer->insert(ta_pointer->begin()+curr_pos, '1');
				}
				else{
					ta_pointer->insert(ta_pointer->begin()+curr_pos, '0');
				}
				curr_pos+=temp_x_size+1;
				temp_last_position++;
			}
			temp_x_size++;
		}
		void remove_site(int r_position){
			int r_curr_pos=r_position;
			for(int i_rm_seq=0; i_rm_seq<y_size; i_rm_seq++){
				ta_pointer->erase(ta_pointer->begin()+r_curr_pos);
				temp_last_position--;
				r_curr_pos+=temp_x_size-1;
			}
			temp_x_size--;
		}
};

class population{
	//
	vector <int> males;
	int n_males=0;
	vector <int> females;
	int n_females=0;
	vector <int> temp_males;
	int n_temp_males=0;
	vector <int> temp_females;
	int n_temp_females=0;
	vector <int> *m_pointer;
	vector <int> *f_pointer;
	vector <int> *tm_pointer;
	vector <int> *tf_pointer;
	public:
	population(vector <char> &initial_sexes){
		for(int i_sex=0; i_sex<initial_sexes.size(); i_sex++){
			if(initial_sexes[i_sex]=='m'){
				males.emplace_back(i_sex*2);
				n_males++;
			}
			else{
				females.emplace_back(i_sex*2);
				n_females++;
			}
		}
	//	cout<<"Initial males:"<<n_males<<endl;
	//	cout<<"Initial females"<<n_females<<endl;
		m_pointer=&males;
		f_pointer=&females;
		tm_pointer=&temp_males;
		tf_pointer=&temp_females;
		//reserve memory for vectors...
		males.reserve(initial_sexes.size());
		females.reserve(initial_sexes.size());
		temp_females.reserve(initial_sexes.size());
		temp_males.reserve(initial_sexes.size());
		/*
		cout<<"Pop n males: "<<n_males<<", male vector:"<<endl;
		for(int i_m=0; i_m<males.size(); i_m++){
			cout<<males[i_m]<<";";
		}
		cout<<endl;
		cout<<"Pop n females: "<<n_females<<", female vector:"<<endl;
		for(int i_f=0; i_f<females.size(); i_f++){
			cout<<females[i_f]<<";";
		}
		cout<<endl;
		*/
	}


		/*
	~population(){
		m_pointer=nullptr;
		f_pointer=nullptr;
		tm_pointer=nullptr;
		tf_pointer=nullptr;

	}
*/

	//rec mode: 0 -> total rec probs at all sites, 1 -> rec probs in blocks between variable sites, 3 -> random

	void make_offspring_gen(genotype_matrix &gen_matrix, const int &og_window_size, const double &og_min_ident, mt19937 &o_randgen, int rec_mode, int o_size, double &og_m_rec_prob, double &og_f_rec_prob ){
		//cout<<"Making offspring gen.."<<endl;
		//cout<<"Making offspring gen, n_males="<<n_males<<",n_females="<<n_females<<endl;
		//iterate over individuals in regular ind array...
		for(int i_o_ind=0; i_o_ind<o_size; i_o_ind++){
			//where in temp array will be their sequences
			int i_o_array_pos=i_o_ind*2;
			//draw mother
		//	cout<<"Adding maternal gamete..."<<endl;
		//	cout<<"Drawing parents..."<<endl;
			int o_mother=draw_unii(0, n_females-1, o_randgen);
			int o_father=draw_unii(0, n_males-1, o_randgen);
		//	cout<<"Done drawing parents..."<<endl;
		//	cout<<"Drew mother, index "<<o_mother<<", positions in genotype array: "<<(*f_pointer)[o_mother]<<","<<(*f_pointer)[o_mother]+1<<endl;
			//TODO: draw recombination position
			
			//int ind_mo_rec_pos=draw_unii(0, 99999);

			//add gamete
			pair <int, int> f_gamete;
			if(draw_unii(0,1, o_randgen)==0){
				f_gamete=make_pair((*f_pointer)[o_mother], (*f_pointer)[o_mother]+1);
			}
			else{
				f_gamete=make_pair((*f_pointer)[o_mother]+1, (*f_pointer)[o_mother]);

			}
			pair <int, int> m_gamete;
			if(draw_unii(0,1, o_randgen)==0){
				m_gamete=make_pair((*m_pointer)[o_father], (*m_pointer)[o_father]+1);
			}
			else{
				m_gamete=make_pair((*m_pointer)[o_father]+1, (*m_pointer)[o_father]);

			}

			if(rec_mode==0){
				gen_matrix.add_recombined_gamete(f_gamete, og_window_size, og_min_ident, o_randgen);
				gen_matrix.add_recombined_gamete(m_gamete, og_window_size, og_min_ident, o_randgen);
			}
			else if(rec_mode==1){
				gen_matrix.add_recombined_gamete_variable_sites(f_gamete, og_window_size ,og_min_ident, o_randgen, og_f_rec_prob);	
				gen_matrix.add_recombined_gamete_variable_sites(m_gamete, og_window_size, og_min_ident, o_randgen, og_m_rec_prob);	
			}
			else if(rec_mode==2){
				gen_matrix.add_recombined_gamete_random(f_gamete, o_randgen);	
				gen_matrix.add_recombined_gamete_random(m_gamete, o_randgen);	
			}
		//	gen_matrix.print_temp_array();
		//	gen_matrix.print_position_array();
			//draw father
		//	cout<<"Adding paternal gamete, n_males="<<n_males<<endl;
		//	cout<<"Drew father, index "<<o_father<<", positions in genotype array: "<<(*m_pointer)[o_father]<<","<<(*m_pointer)[o_father]+1<<endl;
			//TODO: draw recombination position
			//int ind_fa_rec_pos=draw_unii(0, 99999);
			//add gamete
		//	gen_matrix.print_temp_array();
		//	gen_matrix.print_position_array();
			//determine sex
			char o_sex=gen_matrix.determine_sex(make_pair(i_o_array_pos, i_o_array_pos+1));
			//add sequences to appropriate temp array
			if (o_sex=='m'){
				if(tm_pointer->size()<=n_temp_males || tm_pointer->size()==0){
					tm_pointer->emplace_back(i_o_array_pos);
				//	cout<<"Emplacing back male..."<<endl;
				}
				else{
					(*tm_pointer)[n_temp_males]=i_o_array_pos;
				//	cout<<"Overwriting male..."<<endl;
				}
				n_temp_males++;
		//		cout<<"male added to temp array..."<<endl;
			}
			else if (o_sex=='f'){
				if(tf_pointer->size()<=n_temp_females || tf_pointer->size()==0){
					tf_pointer->emplace_back(i_o_array_pos);
				//	cout<<"Emplacing back female..."<<endl;
				}
				else{
					(*tf_pointer)[n_temp_females]=i_o_array_pos;
				//	cout<<"Overwriting female..."<<endl;
				}
				n_temp_females++;
		//		cout<<"female added to temp array..."<<endl;
			}
			else{
				cout<<"Oh dear..."<<endl;
			}
			//gen_matrix.print_temp_array();
		}
		//cut temp array to appropriate size
		//-> works...
		//cout<<"Temp Gen males: "<<n_temp_males<<" : "<<tm_pointer->size()<<endl;
		//cout<<"Temp Gen females: "<<n_temp_females<<" : "<<tf_pointer->size()<<endl;
		while(tm_pointer->size()>n_temp_males){
			tm_pointer->pop_back();
		}
		while(tf_pointer->size()>n_temp_females){
			tf_pointer->pop_back();
		}	
		
		//cout<<"Done adding offspring to temp array..."<<endl;
		//also swap m/f temp and regular arrays...
		swap(m_pointer, tm_pointer);
		n_males=n_temp_males;
		n_temp_males=0;
		swap(f_pointer, tf_pointer);
		n_females=n_temp_females;
		n_temp_females=0;
		//cout<<"Gen males: "<<n_males<<" : "<<m_pointer->size()<<endl;
		//cout<<"Gen females: "<<n_females<<" : "<<f_pointer->size()<<endl;
		//cout<<"Done making offspring, new population has "<<n_males<<" males and "<<n_females<<" females."<<endl;
		//cout<<"Done making offspring gen.."<<endl;
	}
	vector <int> *return_male_pointer(){
		return m_pointer;
	}
	vector <int> *return_female_pointer(){
		return f_pointer;
	}
	vector <int> *return_temp_male_pointer(){
		return tm_pointer;
	}
	vector <int> *return_temp_female_pointer(){
		return tf_pointer;
	}

};

void read_config_file(string cfg_file_name, int &c_ngens, int &c_ninds, int &c_gen_size, double &c_mut_rate, double &c_mr_male_bias, int &c_sd_pos, int &c_rc_wind, double &c_min_ident, int &c_sstat_wind, int &c_sstat_gen, string &c_fst_ofile, string &c_dxy_ofile, string &c_dxx_ofile, double &c_m_rec_prob, double &c_f_rec_prob){
	ifstream config_file;
	config_file.open(cfg_file_name);
	string gc_line;
	while(getline(config_file, gc_line)){
		if(gc_line.length()>0 && gc_line[0]!='#'){
			stringstream config_option(gc_line);
			string param;
			config_option>>param;
			if(param=="ngens"){
				config_option>>c_ngens;
			}
			else if(param=="pop_size"){
				config_option>>c_ninds;
			}
			else if(param=="genome_size"){
				config_option>>c_gen_size;
			}
			else if(param=="mut_rate"){
				config_option>>c_mut_rate;
			}
			else if(param=="mut_mb"){
				config_option>>c_mr_male_bias;
			}
			else if(param=="m_rec_prob"){
				config_option>>c_m_rec_prob;
			}
			else if(param=="f_rec_prob"){
				config_option>>c_f_rec_prob;
			}
			else if(param=="sd_pos"){
				config_option>>c_sd_pos;
			}
			else if(param=="rc_wind"){
				config_option>>c_rc_wind;
			}
			else if(param=="min_ident"){
				config_option>>c_min_ident;
			}
			else if(param=="sstat_wind"){
				config_option>>c_sstat_wind;
			}
			else if(param=="sstat_gen"){
				config_option>>c_sstat_gen;
			}
			else if(param=="fst_ofile"){
				config_option>>c_fst_ofile;
			}
			else if(param=="dxy_ofile"){
				config_option>>c_dxy_ofile;
			}
			else if(param=="dxx_ofile"){
				config_option>>c_dxx_ofile;
			}
		}
	}
}



int main(int argc, char *argv[]){
	//initial parameters:
	//number of generations
	int ngens=1000;
	//population size
	int n_inds=100;
	//window size for effect of divergence in meiosis
	int wind_size=10000;
	//minimum identity to recombination prob of block to 0
	double min_ident=0.6;
	//genome size
	int gen_size=100000;
	//mutation rate
	double mutation_rate=2.5e-7;
	//male bias of mutation rate
	double mr_male_bias=1;
	//probability of a single crossing over in male meiosis
	double m_rec_prob=1;
	//probability of a single crossing over in female meiosis
	double f_rec_prob=1;
	//sd position
	int sd_loc=50000;
	//mean window size for summary statistics
	int mean_rp_wind_size=5000;
	//how often shall we print summary stats?(every? generations)
	int sum_stat=10;
	//sumstats m/f FST output
	string fst_ofile="m_f_fst.out";
	//dxx sumstat file
	string dxx_ofile="dxx.out";
	//dxy sumstat file	
	string dxy_ofile="dxy.out";


	//read config file
	read_config_file(argv[1], ngens, n_inds, gen_size, mutation_rate, mr_male_bias, sd_loc, wind_size, min_ident, mean_rp_wind_size, sum_stat, fst_ofile, dxy_ofile, dxx_ofile, m_rec_prob, f_rec_prob );

	//open FST sumstat output
	ofstream fst_out;
	ofstream dxx_out;
	ofstream dxy_out;
	dxx_out.open(dxx_ofile);
	dxy_out.open(dxy_ofile);
	fst_out.open(fst_ofile);
	fst_out<<"generation\tN_variants";
	dxx_out<<"generation";
	dxy_out<<"generation";
	for(int i_fst_h=mean_rp_wind_size; i_fst_h<=gen_size; i_fst_h+=mean_rp_wind_size){
		fst_out<<"\t"<<i_fst_h;
		dxx_out<<"\t"<<i_fst_h;
		dxy_out<<"\t"<<i_fst_h;
	}
	fst_out<<endl;
	dxx_out<<endl;
	dxy_out<<endl;
	random_device rand_dev;
        mt19937 r_generator(rand_dev());
	vector <char> initial_sex_vector;
	initial_sex_vector.reserve(n_inds);
	//assign intial females
	initial_sex_vector.assign(n_inds/2, 'f');
	//add males
	for(int i_m=0; i_m<n_inds/2; i_m++){
		initial_sex_vector.emplace_back('m');
	}
	genotype_matrix initial_matrix(initial_sex_vector, gen_size, sd_loc, mean_rp_wind_size );
	population initial_pop(initial_sex_vector);
	for(int i_gen=0; i_gen<ngens; i_gen++){
		initial_pop.make_offspring_gen(initial_matrix, wind_size, min_ident, r_generator, 1, n_inds, m_rec_prob, f_rec_prob );
//		cout<<"Made offspring..."<<endl;
		initial_matrix.temp_remove_fixed_sites();
		initial_matrix.temp_add_mutations(mutation_rate, mr_male_bias, initial_pop.return_temp_female_pointer(), initial_pop.return_temp_male_pointer(), r_generator);
		initial_matrix.swap_arrays();
//		initial_matrix.print_array();
		if(i_gen%sum_stat==0){
			cout<<"Generation "<<i_gen<<endl;
			initial_matrix.write_window_fst(mean_rp_wind_size, initial_pop.return_male_pointer(), initial_pop.return_female_pointer(), i_gen, fst_out);
			initial_matrix.write_window_dxy(mean_rp_wind_size, initial_pop.return_female_pointer(), initial_pop.return_male_pointer(), dxx_out, dxy_out, i_gen);
		//	initial_matrix.print_array();
		}

	}
	fst_out.close();
	dxx_out.close();
	dxy_out.close();
return 0;
};
