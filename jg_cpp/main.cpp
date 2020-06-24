#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <stdlib.h>
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
	//maybe just initialize once!
	//random_device w_rand_dev;
        //mt19937 w_generator(w_rand_dev());
	//until here
	discrete_distribution<int> w_distribution (w_prob_vector.begin(), w_prob_vector.end());
	//cout<<"Done making weighted draw"<<endl;
	return w_distribution(w_generator);
}
int draw_poisson(double &p_mean, mt19937 &p_generator){
	poisson_distribution<int> p_distribution (p_mean);
	return p_distribution(p_generator);
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
	
	
	//mean rec_prob vector:
	vector <double> mean_rec_prob_vector;
	//takes vector of sexes of individuals and adds them in this order to the genotype array...
	//determines left positions in for block var array and recombination probs
	vector <int> recomb_left_positions_var;		
	vector <double> recomb_probs_var;
	vector <int> block_var;	
	//probability vectors for hotspot model
	vector <int> hs_rec_sites;
	vector <double> hs_original_rec_probs_m;
	vector <double> hs_original_rec_probs_f;

	vector <double> hs_modified_rec_probs;
	vector <int> hs_left_rec_sites;

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
		genotype_array.reserve(100000);
		temp_array.reserve(100000);
		ga_pointer=&genotype_array;
		ta_pointer=&temp_array;
		prob_vector.reserve(gen_size-1);
		prob_vector.assign(gen_size-1,1);
		//number of windows in mean rec prob vector
		//int mean_rec_prob_vector_size=gen_size/mean_rec_prob_window_size;
		//mean_rec_prob_vector.reserve(mean_rec_prob_vector_size);
		//mean_rec_prob_vector.assign(mean_rec_prob_vector_size,0);
		recomb_left_positions_var.reserve(10000);
		recomb_probs_var.reserve(10000);
	}
	~genotype_matrix(){
		ga_pointer=nullptr;
		ta_pointer=nullptr;
	}
	void read_hotspot_conf(string &hs_conf_file){
		cout<<"Reading hotspot conf"<<endl;
		ifstream hs_config_file;
		hs_config_file.open(hs_conf_file);
		string hs_line;
		int hs_position;
		double hs_male_prob;
		double hs_female_prob;
		while(getline(hs_config_file, hs_line)){
			if(hs_line.length()>0 && hs_line[0]!='#'){
				stringstream hs_config_option(hs_line);
				hs_config_option>>hs_position;
				if(hs_position<=genome_size){
					hs_rec_sites.emplace_back(hs_position);
					hs_config_option>>hs_male_prob;
					hs_original_rec_probs_m.emplace_back(hs_male_prob);
					hs_config_option>>hs_female_prob;
					hs_original_rec_probs_f.emplace_back(hs_female_prob);
					hs_modified_rec_probs.emplace_back(hs_female_prob);
					hs_left_rec_sites.emplace_back(0);
				}
				else{
					cout<<"Warning: recombination site "<<hs_position<<" is larger than the genome size. Ignoring site."<<endl;
				}
			}
		}
		//add prob for unsimulated region
		hs_modified_rec_probs.emplace_back(0);
		cout<<"Done reading hotspot.conf. Hotspot vector is "<<hs_original_rec_probs_m.size()<<" long."<<endl;
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
	
	double determine_block_rec_prob(const int &brp_block_size, const int &brp_block_cov, const int &brp_wind_size, const double &brp_min_identity, const double &brp_var_effect){
		double brp_rec_prob=0;
		if(brp_block_size>0 && double(brp_wind_size-brp_block_cov)/brp_wind_size>brp_min_identity && brp_var_effect*brp_block_cov<brp_wind_size){
			//brp_rec_prob=(brp_block_size*(double(brp_wind_size-(brp_var_effect*brp_block_cov))/brp_wind_size))/genome_size;
			brp_rec_prob=(brp_block_size*(1-(double(brp_var_effect*brp_block_cov)/brp_wind_size)))/genome_size;
			//cout<<"Making recombination probability for dynamic model. Prob: "<<brp_rec_prob<<", var_effect: "<<brp_var_effect<<" ,n_variants: "<<brp_block_cov<<", block size: "<<brp_block_size<<", window_size: "<<brp_wind_size<<", genome_size: "<<genome_size<<endl;
			//TODO: continue here...test var_effect...
		}
		return brp_rec_prob;
	}
	

	void recomb_block(const int &rb_window_size, int &rb_curr_sum_var, const int &rb_var_site_count, int &rb_prev_block, int &rb_curr_block, const int &rb_curr_block_cov, int &rb_n_blocks, double &rb_curr_block_sum,const double &rb_min_identity, const double &rb_var_effect ){
		double curr_block_prob;
//		cout<<"Calling recomb_block, curr_block="<<rb_curr_block<<", prev_block="<<rb_prev_block<<endl;
		if(rb_curr_sum_var<rb_var_site_count){
			while(rb_curr_sum_var<rb_var_site_count-1 && rb_prev_block>position_array[recomb_left_positions_var[rb_curr_sum_var]]){
				/*
				if(recomb_probs_var.size()>rb_curr_sum_var){
					recomb_probs_var[rb_curr_sum_var]=rb_curr_block_sum;
				}
				else{
					recomb_probs_var.emplace_back(rb_curr_block_sum);
				}
				*/
				recomb_probs_var.emplace_back(rb_curr_block_sum);
				rb_curr_block_sum=0;
				if(rb_curr_block_sum<rb_var_site_count-1){
					rb_curr_sum_var++;
				}
			}
			if(rb_curr_block<=position_array[recomb_left_positions_var[rb_curr_sum_var]]){
//				cout<<"\t1 Adding block prob. Begin: "<<rb_prev_block<<", End: "<<rb_curr_block<<", Cov: "<<rb_curr_block_cov<<endl;
		/*		
				//only add sum if block identity is larger than min identity...
				if((double(rb_curr_block_cov)/double(rb_curr_block-rb_prev_block))<1-rb_min_identity){
				//include min coverage and site effect here...
					curr_block_prob=(rb_curr_block-rb_prev_block-(rb_var_effect*(rb_curr_block_cov/double(rb_curr_block-rb_prev_block))))/genome_size;
					if(curr_block_prob<0){
						curr_block_prob=0;
					}
					rb_curr_block_sum+=curr_block_prob;
				}
			*/
				rb_curr_block_sum+=determine_block_rec_prob(rb_curr_block-rb_prev_block, rb_curr_block_cov,rb_window_size, rb_min_identity, rb_var_effect );
			}
				
			else{
//				cout<<"\t2 Adding block prob. Begin: "<<rb_prev_block<<", End: "<<position_array[recomb_left_positions_var[rb_curr_sum_var]]<<", Cov: "<<rb_curr_block_cov<<endl;
				/*
				if(position_array[recomb_left_positions_var[rb_curr_sum_var]]>rb_prev_block && (double(rb_curr_block_cov)/double(position_array[recomb_left_positions_var[rb_curr_sum_var]]-rb_prev_block))<1-rb_min_identity){
					curr_block_prob=(position_array[recomb_left_positions_var[rb_curr_sum_var]]-rb_prev_block-(rb_var_effect*(rb_curr_block_cov/double(position_array[recomb_left_positions_var[rb_curr_sum_var]]-rb_prev_block))))/genome_size;
					if(curr_block_prob<0){
						curr_block_prob=0;
					}
					rb_curr_block_sum+=curr_block_prob;
				}
				*/
				rb_curr_block_sum+=determine_block_rec_prob(position_array[recomb_left_positions_var[rb_curr_sum_var]]-rb_prev_block, rb_curr_block_cov,rb_window_size, rb_min_identity, rb_var_effect );


				recomb_probs_var.emplace_back(rb_curr_block_sum);
				
				//continue here...
				while(rb_curr_sum_var+1<rb_var_site_count && position_array[recomb_left_positions_var[rb_curr_sum_var+1]]<rb_curr_block){
//					cout<<"\t3 Adding block prob. Begin: "<<position_array[recomb_left_positions_var[rb_curr_sum_var]]+1<<", End: "<<position_array[recomb_left_positions_var[rb_curr_sum_var+1]]<<", Cov: "<<rb_curr_block_cov<<endl;
					/*
					if((double(rb_curr_block_cov)/double(position_array[recomb_left_positions_var[rb_curr_sum_var+1]]-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1))<1-rb_min_identity){
					
						rb_curr_block_sum=(position_array[recomb_left_positions_var[rb_curr_sum_var+1]]-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1-(rb_var_effect*(rb_curr_block_cov/double(position_array[recomb_left_positions_var[rb_curr_sum_var+1]]-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1))))/genome_size;
						if(rb_curr_block_sum<0){
							rb_curr_block_sum=0;
						}
					}
					else{
						rb_curr_block_sum=0;
					}
					*/
					rb_curr_block_sum=determine_block_rec_prob(position_array[recomb_left_positions_var[rb_curr_sum_var+1]]-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1, rb_curr_block_cov,rb_window_size, rb_min_identity, rb_var_effect );


					
					/*
					if(recomb_probs_var.size()>rb_curr_sum_var){
						recomb_probs_var[rb_curr_sum_var]=rb_curr_block_sum;
					}
					else{
						recomb_probs_var.emplace_back(rb_curr_block_sum);
					}
					*/
					recomb_probs_var.emplace_back(rb_curr_block_sum);
					if(rb_curr_sum_var<rb_var_site_count){
						rb_curr_sum_var++;
					}
				}
//				cout<<"\t4 Adding block prob. Begin: "<<position_array[recomb_left_positions_var[rb_curr_sum_var]]+1<<", End: "<<rb_curr_block<<", Cov: "<<rb_curr_block_cov<<endl;
				/*
				if((double(rb_curr_block_cov)/double(rb_curr_block-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1))<1-rb_min_identity){
					rb_curr_block_sum=(rb_curr_block-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1-(rb_var_effect*(rb_curr_block_cov/double(rb_curr_block-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1))))/genome_size;
					if(rb_curr_block_sum<0){
						rb_curr_block_sum=0;
					}
				}
				else{
					rb_curr_block_sum=0;
				}
				*/
				rb_curr_block_sum=determine_block_rec_prob(rb_curr_block-position_array[recomb_left_positions_var[rb_curr_sum_var]]+1, rb_curr_block_cov,rb_window_size, rb_min_identity, rb_var_effect );

				if(rb_curr_sum_var<rb_var_site_count){
					rb_curr_sum_var++;
				}
			}
		}
		else{
//			cout<<"\t5 Adding block prob. Begin: "<<rb_prev_block<<", End: "<<rb_curr_block<<", Cov: "<<rb_curr_block_cov<<endl;
			/*
			if(double(rb_curr_block_cov)/double(rb_curr_block-rb_prev_block)<1-rb_min_identity){
				curr_block_prob=(rb_curr_block-rb_prev_block-(rb_var_effect*(rb_curr_block_cov/double(rb_curr_block-rb_prev_block))))/genome_size;
				if(curr_block_prob<0){
					curr_block_prob=0;
				}
				rb_curr_block_sum+=curr_block_prob;
			}
			
		
		*/
			rb_curr_block_sum+=determine_block_rec_prob(rb_curr_block-rb_prev_block, rb_curr_block_cov,rb_window_size, rb_min_identity, rb_var_effect );
		}	
		rb_prev_block=rb_curr_block+1;
		rb_n_blocks++;
	}


	void add_recombined_gamete_variable_sites(pair <int, int> ind_array_pos_var, const int &rg_window_size_var, const double &rg_min_ident_var, const double &rg_var_effect, mt19937 &rc_randgen_var, const double &genome_rec_prob, int &rg_rec_point,vector <double> **rec_vector_probs_var,  vector <int> **rec_vector_positions_var,  vector <int> **rec_vector_pos_array){
			//vector <double> recomb_probs_var;
			//vector <int> recomb_left_positions_var;
			int window_end_it=0;
			int left_window;
			int curr_block_cov=0;
			int n_var_sites=0;
			int n_blocks=0;
			int curr_block=0;
			int prev_block=0;
			double curr_block_sum=0;
			int curr_sum_var=0;
			double add_block;
			//empty recomb_probs_var
			recomb_probs_var.clear();
			//iterate over position array
			for(int i_rc_pos=0; i_rc_pos<position_array.size(); i_rc_pos++){
				//only use sides that are variable between both seqs
				if((*ga_pointer)[(ind_array_pos_var.first*x_size)+i_rc_pos] != (*ga_pointer)[(ind_array_pos_var.second*x_size)+i_rc_pos]){
					//remember variable site first!
					if(recomb_left_positions_var.size()>n_var_sites){
						recomb_left_positions_var[n_var_sites]=i_rc_pos;
					}
					else{
						recomb_left_positions_var.emplace_back(i_rc_pos);
					}
					left_window=position_array[i_rc_pos]-(rg_window_size_var/2);	
					while(position_array[recomb_left_positions_var[window_end_it]]+(rg_window_size_var/2)<left_window){
						curr_block=position_array[recomb_left_positions_var[window_end_it]]+(rg_window_size_var/2);
						if(prev_block<curr_block){
//							cout<<"rc_call 1"<<endl;
							recomb_block(rg_window_size_var, curr_sum_var, n_var_sites, prev_block, curr_block, curr_block_cov, n_blocks, curr_block_sum, rg_min_ident_var, rg_var_effect );
						}
						window_end_it++;
						curr_block_cov--;
					}
					if(left_window>0){
						if(n_blocks>0){
							if(prev_block<left_window){
//								cout<<"rc_call 2"<<endl;
								recomb_block(rg_window_size_var,curr_sum_var, n_var_sites, prev_block, left_window, curr_block_cov, n_blocks, curr_block_sum, rg_min_ident_var, rg_var_effect );
							}

						}
						else{
//							cout<<"rc_call 3"<<endl;
							recomb_block(rg_window_size_var, curr_sum_var, n_var_sites, prev_block, left_window, curr_block_cov, n_blocks, curr_block_sum, rg_min_ident_var, rg_var_effect );
						}
					}
					curr_block_cov++;
					n_var_sites++;
				}
			}
	
			//test if both sequences are identical
			if(n_var_sites>0){
				//finish looping window_end_it
				while(window_end_it<n_var_sites && position_array[recomb_left_positions_var[window_end_it]]+(rg_window_size_var/2)<genome_size){
					curr_block=position_array[recomb_left_positions_var[window_end_it]]+(rg_window_size_var/2);
					//continue here			
					if(prev_block<curr_block){
//						cout<<"rc_call 4"<<endl;
						recomb_block(rg_window_size_var, curr_sum_var, n_var_sites, prev_block, curr_block, curr_block_cov, n_blocks, curr_block_sum, rg_min_ident_var, rg_var_effect );
					}
					window_end_it++;
					curr_block_cov--;
				}
				//finish looping over var_sites
				int prev_position=prev_block;
				while(curr_sum_var<n_var_sites){
//					cout<<"\t6 Adding block prob. Begin: "<<prev_position+1<<", End: "<<position_array[recomb_left_positions_var[curr_sum_var]]<<", Cov: "<<curr_block_cov<<endl;
						/*
					if(position_array[recomb_left_positions_var[curr_sum_var]]>prev_position && (double(curr_block_cov)/double(position_array[recomb_left_positions_var[curr_sum_var]]-prev_block))<1-rg_min_ident_var){
						add_block=(position_array[recomb_left_positions_var[curr_sum_var]]-prev_position-(rg_var_effect*(curr_block_cov/double(position_array[recomb_left_positions_var[curr_sum_var]]-prev_position))))/genome_size;
						if(add_block<0){
							add_block=0;
						}
						curr_block_sum+=add_block;
					}
						*/
						curr_block_sum+=determine_block_rec_prob(position_array[recomb_left_positions_var[curr_sum_var]]-prev_position+1  , curr_block_cov, rg_window_size_var, rg_min_ident_var, rg_var_effect );

					/*
					if(recomb_probs_var.size()>curr_sum_var){
						recomb_probs_var[curr_sum_var]=curr_block_sum;
					}
					else{
						recomb_probs_var.emplace_back(curr_block_sum);
					}
					*/
					recomb_probs_var.emplace_back(curr_block_sum);
					prev_position=position_array[recomb_left_positions_var[curr_sum_var]];
					curr_block_sum=0;
					curr_sum_var++;
				}
		
				//Add final Block...
				if(prev_block<genome_size-1){
					if(prev_block<position_array[recomb_left_positions_var[n_var_sites-1]]){
//						cout<<"\t7 Adding block prob. Begin: "<<position_array[recomb_left_positions_var[n_var_sites-1]]+1<<", End: "<<genome_size-1<<", Cov: "<<curr_block_cov<<endl;
						/*
						if((double(curr_block_cov)/double(genome_size-1-position_array[recomb_left_positions_var[n_var_sites-1]]))<1-rg_min_ident_var){
							add_block=(genome_size-1-position_array[recomb_left_positions_var[n_var_sites-1]]-(rg_var_effect*(curr_block_cov/double(genome_size-1-position_array[recomb_left_positions_var[n_var_sites-1]]))))/genome_size;
							if(add_block<0){
								add_block=0;
							}
							curr_block_sum+=add_block;
						}
							*/
							curr_block_sum+=determine_block_rec_prob(genome_size-position_array[recomb_left_positions_var[n_var_sites-1]], curr_block_cov, rg_window_size_var, rg_min_ident_var, rg_var_effect );

					}
					else{
//						cout<<"\t8 Adding block prob. Begin: "<<prev_block<<", End: "<<genome_size-1<<", Cov: "<<curr_block_cov<<endl;
						/*
						if((double(curr_block_cov)/double(genome_size-1-prev_block))<1-rg_min_ident_var){
							add_block=(genome_size-1-prev_block-(rg_var_effect*double(curr_block_cov/(genome_size-1-prev_block))))/genome_size;
							if(add_block<0){
								add_block=0;
							}
							curr_block_sum+=add_block;
						}
						*/
						curr_block_sum+=determine_block_rec_prob(genome_size-1-prev_block, curr_block_cov, rg_window_size_var, rg_min_ident_var, rg_var_effect );

					}
					/*
					if(recomb_probs_var.size()>curr_sum_var){
						recomb_probs_var[curr_sum_var]=curr_block_sum;
					}
					else{
						recomb_probs_var.emplace_back(curr_block_sum);
					}
					*/
					recomb_probs_var.emplace_back(curr_block_sum);
			
				}
				/*
				//cut back prob vector
				while(recomb_probs_var.size()>n_var_sites+1){
					recomb_probs_var.pop_back();
				}
				*/
				//add genomic probability to last window
				//recomb_probs_var.back()+=genome_rec_prob;
				recomb_probs_var.emplace_back(genome_rec_prob);
				/*
				double total_prob=0;
				for(int i_rpv=0; i_rpv<recomb_probs_var.size(); i_rpv++){
					cout<<recomb_probs_var[i_rpv]<<" ";
					total_prob+=recomb_probs_var[i_rpv];
				}
				cout<<", Prob sum: "<<total_prob<<endl;
				assert(total_prob<=genome_rec_prob+1);
				//make weighted draw
				cout<<"Making weighted draw. Prob vector size: "<<recomb_probs_var.size()<<",N var sites: "<<n_var_sites<<", recomb_left_positions_var size: "<<recomb_left_positions_var.size()<<endl;
				*/
				int rc_p_index=make_weighted_draw(recomb_probs_var, rc_randgen_var);
				


				//check if position is outside the simulated region or at the beginning or at the end 
				if(rc_p_index>=recomb_probs_var.size()-2 || rc_p_index==0 ){
					//make non recombinant gamete
					if(draw_unii(0,1, rc_randgen_var)==0){
						non_recombined_gamete(ind_array_pos_var.first);
					}
					else{
						non_recombined_gamete(ind_array_pos_var.second);
					}
				}
				//otherwise make recombinant gametes...
				else{	
					//I think this was wrong!			
					//recombined_gamete(recomb_left_positions_var[rc_p_index-1], ind_array_pos_var);
					recombined_gamete(recomb_left_positions_var[rc_p_index-1], ind_array_pos_var);
				}
				//
				if(rc_p_index==recomb_probs_var.size()-1){
					rg_rec_point=genome_size+1;
				}
				/*
				else if(){
				
				
				}
				*/
				else{
				
					rg_rec_point=position_array[recomb_left_positions_var[rc_p_index]];
				}
			}
			else{
				//make non-recombinant gamete
				non_recombined_gamete(ind_array_pos_var.first);
			}
			*rec_vector_positions_var=&recomb_left_positions_var;
			*rec_vector_probs_var=&recomb_probs_var;
			*rec_vector_pos_array=&position_array;
			
	}
	void add_recombined_gamete_hotspot(pair <int, int> ind_array_pos_hs, const int &hs_window_size, const double &hs_min_ident, const double &hs_var_effect, mt19937 &hs_randgen_var, const double &hs_genome_rec_prob, bool hs_male_rec, int &hs_rec_point, vector <double> **rec_vector_probs_hs, vector <int> **rec_vector_positions_hs){
		int last_rec_pos=0;
		int curr_rec_pos=0;
		int hs_left_rec_pos=0;
		int hs_window_var_site_count=0;
		bool found_left_rec_pos=false;
		for(int i_hs_pos=0; i_hs_pos<hs_rec_sites.size(); i_hs_pos++){
			//advance last_rec_prob until it's within window size
			while(position_array[last_rec_pos]<hs_rec_sites[i_hs_pos]-(hs_window_size/2) && position_array[last_rec_pos+1]<hs_rec_sites[i_hs_pos]+(hs_window_size/2)){
				last_rec_pos++;
			}
			//count variable sites in window
			curr_rec_pos=last_rec_pos;
			//
			hs_left_rec_pos=curr_rec_pos;
			while(position_array[curr_rec_pos]<hs_rec_sites[i_hs_pos]+(hs_window_size/2)){
				//remember where rec spots are...
				if(position_array[curr_rec_pos]==hs_rec_sites[i_hs_pos]){
					hs_left_rec_pos=curr_rec_pos;
				//	cout<<"1: Found left_rec_pos, site: "<<hs_rec_sites[i_hs_pos]<<endl;
					found_left_rec_pos=true;		
				} 
				else if(found_left_rec_pos==false && position_array[curr_rec_pos]>hs_rec_sites[i_hs_pos] && curr_rec_pos>0){
					hs_left_rec_pos=curr_rec_pos-1;
				//	cout<<"2: Found left_rec_pos, site: "<<hs_rec_sites[i_hs_pos]<<endl;
					found_left_rec_pos=true;
				}
				if((*ga_pointer)[(ind_array_pos_hs.first*x_size)+curr_rec_pos] != (*ga_pointer)[(ind_array_pos_hs.second*x_size)+curr_rec_pos]){
					hs_window_var_site_count++;
				}
				curr_rec_pos++;
			}
			if(found_left_rec_pos==false && position_array[curr_rec_pos]<hs_rec_sites[i_hs_pos])
			{
			//	cout<<"3: Found left_rec_pos, site: "<<hs_rec_sites[i_hs_pos]<<endl;
				hs_left_rec_pos=curr_rec_pos;
			}
			hs_left_rec_sites[i_hs_pos]=hs_left_rec_pos;
			found_left_rec_pos=false;
			double hs_mod_rec_prob;
			if(hs_male_rec==true){
				if(1-(hs_window_var_site_count/hs_window_size)>hs_min_ident){
					hs_mod_rec_prob=hs_original_rec_probs_m[i_hs_pos]*(1-((hs_var_effect*hs_window_var_site_count)/hs_window_size));
					//cout<<"Making male recombination probability for hotspot model. Prob: "<<hs_mod_rec_prob<<", var_effect: "<<hs_var_effect<<" ,n_variants: "<<hs_window_var_site_count<<", a priori rec prob: "<<hs_original_rec_probs_m[i_hs_pos]<<", window_size: "<<hs_window_size<<", genome_size: "<<genome_size<<endl;
				}
				else{
					hs_mod_rec_prob=0;
				}
			}
			else{
				if(1-(hs_window_var_site_count/hs_window_size)>hs_min_ident){
					hs_mod_rec_prob=hs_original_rec_probs_f[i_hs_pos]*(1-((hs_var_effect*hs_window_var_site_count)/hs_window_size));
					//cout<<"Making female recombination probability for hotspot model. Prob: "<<hs_mod_rec_prob<<", var_effect: "<<hs_var_effect<<" ,n_variants: "<<hs_window_var_site_count<<", a priori rec prob: "<<hs_original_rec_probs_f[i_hs_pos]<<", window_size: "<<hs_window_size<<", genome_size: "<<genome_size<<endl;
				}
				else{
					hs_mod_rec_prob=0;
				}
			}
			if(hs_mod_rec_prob<0){
				hs_mod_rec_prob=0;
			}
			hs_modified_rec_probs[i_hs_pos]=hs_mod_rec_prob;
			hs_window_var_site_count=0;
		}
		//add genome_rec_prob
		hs_modified_rec_probs.back()=hs_genome_rec_prob;
		/*
		print_position_array();
		cout<<"Hs vector:"<<endl;
		for(int i_hs_mod=0; i_hs_mod<hs_modified_rec_probs.size(); i_hs_mod++){
			cout<<"Site: "<<hs_rec_sites[i_hs_mod]<<" Prob: "<<hs_modified_rec_probs[i_hs_mod]<<" Left_rec_sîte: "<<hs_left_rec_sites[i_hs_mod]<<", In position array: "<<position_array[hs_left_rec_sites[i_hs_mod]]<<endl;
		}
		*/
		//make weighted draw
		int hs_p_index=make_weighted_draw(hs_modified_rec_probs, hs_randgen_var);
		if(hs_p_index==hs_modified_rec_probs.size()-1){
			hs_rec_point=genome_size+1;
		}
		else{
			hs_rec_point=hs_rec_sites[hs_p_index];
		}
		//cout<<"hs_index: "<<hs_p_index<<" hs_rec_point: "<<hs_rec_point<<endl;
		//non-recombined gamete
		//if(hs_p_index>=hs_modified_rec_probs.size()-2){
		if(hs_p_index==hs_modified_rec_probs.size()-1){
			if(draw_unii(0,1, hs_randgen_var)==0){
				non_recombined_gamete(ind_array_pos_hs.first);
			}
			else{
				non_recombined_gamete(ind_array_pos_hs.second);
			}
		}
		//recombined gamete
		else{
			recombined_gamete(hs_left_rec_sites[hs_p_index], ind_array_pos_hs);
			
		}
		*rec_vector_probs_hs=&hs_modified_rec_probs;
		*rec_vector_positions_hs=&hs_rec_sites;
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
		double mean_mutations=(mut_males->size()+mut_females->size())*2*genome_size*mutation_rate;
		int total_number_of_mutations=draw_poisson(mean_mutations, am_randgen);
		//cout<<"N_mutations:"<<total_number_of_mutations<<endl;
		//draw total number of mutations in males (use temp female array)
		int female_muts=total_number_of_mutations*(1/(male_bias+1));
		int male_muts=total_number_of_mutations-female_muts;
		//test if position array is large enough
		if(position_array.size()+total_number_of_mutations>genome_size){
			cout<<"Error: there are not enough free sites in the simulated region to add all mutations. Decrease mutation rate or population size or increase Genome size!"<<endl;
			exit(EXIT_FAILURE);	
		}		
		
	
		//draw total number of mutations in females
		//TODO: continue here
		//cout<<"Total muts: "<<total_number_of_mutations<<",male muts: "<<male_muts<<",female_muts: "<<female_muts<<endl;
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
		sdr_file<<sdr_size<<endl;
	}
/*	
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
*/
	//m/f FST
	void write_window_fst(int fst_window_size, vector <int> *fst_m_pointer, vector <int> *fst_f_pointer, const int &fst_gen, ostream &fst_outfile){
		//make vector of windows
		int n_windows=genome_size/fst_window_size;
		//vector <double> fst_per_window(n_windows, 0);
		//just cout for now, later into an output file...
	//	cout<<"FST: ";
		//iterate over position array
		fst_outfile<<fst_gen<<"\t"<<position_array.size();
		

		int curr_window_end=fst_window_size;
		int curr_window_begin=0;
		int i_fst=0;
		int curr_fst_pos=position_array[0];
		//iterate over windows instead
		for(int i_fst_wind=0; i_fst_wind<n_windows; i_fst_wind++){
			double curr_window_fst=0;
			int curr_window_sites=0;

			while(curr_fst_pos<curr_window_end && curr_fst_pos>=curr_window_begin && i_fst<position_array.size()){
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
					
				curr_window_fst+=site_fst;
				curr_window_sites++;
				i_fst++;
				curr_fst_pos=position_array[i_fst];
			}
			if(curr_window_sites>0){
			//	fst_per_window[fst_window_counter]=curr_window_fst/curr_window_sites;
				fst_outfile<<"\t"<<curr_window_fst/double(curr_window_sites);
			}
			else{
			//	fst_per_window[fst_window_counter]=0;
				fst_outfile<<"\tNA";
			}
			curr_window_begin=curr_window_end;
			curr_window_end+=fst_window_size;
		}
		fst_outfile<<endl;
	}
	void write_window_dxy(int dxy_window_size, vector <int> *dxx_f_pointer, vector <int> *dxy_m_pointer, ostream &dxx_outfile, ostream &dxy_outfile, int dxy_gen){
		int curr_window_end=dxy_window_size;
		int curr_window_begin=0;
		int n_windows=genome_size/dxy_window_size;
		dxy_outfile<<dxy_gen;
		dxx_outfile<<dxy_gen;
		int dxy_pos_array_count=0;
		//iterate over windows
		for(int i_dxy_window=0; i_dxy_window<n_windows; i_dxy_window++){
			double curr_window_dxy=0;
			double curr_window_dxx=0;
			while(position_array[dxy_pos_array_count]<curr_window_end && position_array[dxy_pos_array_count]>=curr_window_begin && dxy_pos_array_count<position_array.size()){
				/*
				int local_window_size;

				if(curr_window_end<=genome_size){
					local_window_size=dxy_window_size;
				}
				else{
					local_window_size=curr_window_end-genome_size;
				}
				*/

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
			curr_window_begin=curr_window_end;
			curr_window_end+=dxy_window_size;
		}
		dxy_outfile<<endl;
		dxx_outfile<<endl;
	}
	//, ofstream &rv_ofile 
	/*
	void write_rec_vectors(vector <int> &rv_var_sites, vector <double> &rv_rec_probs){
		cout<<"Rec probs, size="<<rv_rec_probs.size()<<endl;
		for(int i_r_rv=0; i_r_rv<rv_rec_probs.size(); i_r_rv++){
			cout<<rv_rec_probs[i_r_rv]<<";";
		}
		cout<<endl;
		if(rv_var_sites.size()>=rv_rec_probs.size()-2){
			cout<<"var_sites, size="<<rv_var_sites.size()<<endl;
			for(int i_r_vs=0; i_r_vs<rv_rec_probs.size()-2; i_r_vs++){
				cout<<position_array[rv_var_sites[i_r_vs]]<<";";
			}
		}
		cout<<genome_size<<endl<<"------------------------------------------------"<<endl;
	}
	*/
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
	//		cout<<"Adding nonrec gamete! ta_vector size: "<<ta_pointer->size()<<", temp_last_position: "<<temp_last_position<<endl;
	
			int seq_start_position=seq_id*x_size;
			while(ta_pointer->size()<temp_last_position+x_size){
				ta_pointer->emplace_back('0');
	//			cout<<"Emplacing zero to temp vector, new temp vector size: "<<ta_pointer->size()<<endl;
			}


			for(int i_non_rec_pos=seq_start_position; i_non_rec_pos<(seq_id+1)*x_size; i_non_rec_pos++){
			//	cout<<"Copying site "<<i_non_rec_pos<<" from genotype array position "<<i_non_rec_pos<<" to temp array position "<<temp_last_position<<endl;
				(*ta_pointer)[temp_last_position]=(*ga_pointer)[i_non_rec_pos];
				temp_last_position++;
			//	cout<<"Temp last position: "<<temp_last_position<<endl;
			}
			temp_x_size=x_size;
	//		cout<<"Added nonrec gamete..."<<endl;
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

	vector <int> m_rec_points;
	vector <int> f_rec_points;
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

	void make_offspring_gen(genotype_matrix &gen_matrix, const int &og_window_size, const double &og_min_ident, mt19937 &o_randgen, int rec_mode, int o_size, double &og_rec_prob, bool record_mf_rec_points, ofstream &m_rec_vector_out,bool record_m_rec_vector,const double &og_var_effect, int o_genome_size ){
		//empty record_probs
		if(record_mf_rec_points==true){
			m_rec_points.clear();
			f_rec_points.clear();
		}
		vector <double> *m_rec_vector=nullptr;
		vector <int> *m_rec_pos=nullptr;
		vector <int> *m_rec_pos_array=nullptr;

		


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
			int male_rec_point=0;
			int female_rec_point=0;

			if(rec_mode==0){
				gen_matrix.add_recombined_gamete(f_gamete, og_window_size, og_min_ident, o_randgen);
				gen_matrix.add_recombined_gamete(m_gamete, og_window_size, og_min_ident, o_randgen);
			}
			else if(rec_mode==1){

				gen_matrix.add_recombined_gamete_variable_sites(f_gamete, og_window_size ,og_min_ident, og_var_effect,  o_randgen, og_rec_prob, female_rec_point, &m_rec_vector, &m_rec_pos, &m_rec_pos_array);	
				gen_matrix.add_recombined_gamete_variable_sites(m_gamete, og_window_size, og_min_ident, og_var_effect, o_randgen, og_rec_prob, male_rec_point, &m_rec_vector, &m_rec_pos, &m_rec_pos_array);	
			
				if(record_mf_rec_points==true){
					m_rec_points.push_back(male_rec_point);
					f_rec_points.push_back(female_rec_point);
				}
				if(record_m_rec_vector==true){
					write_rec_vector(m_rec_vector, m_rec_pos, m_rec_pos_array, m_rec_vector_out, o_genome_size);	
				}
						
			}
			else if(rec_mode==2){
				gen_matrix.add_recombined_gamete_random(f_gamete, o_randgen);	
				gen_matrix.add_recombined_gamete_random(m_gamete, o_randgen);	
			}

			else if(rec_mode==3){
				gen_matrix.add_recombined_gamete_hotspot(f_gamete, og_window_size, og_min_ident, og_var_effect, o_randgen, og_rec_prob, false, female_rec_point, &m_rec_vector, &m_rec_pos);
				gen_matrix.add_recombined_gamete_hotspot(m_gamete, og_window_size, og_min_ident, og_var_effect, o_randgen, og_rec_prob, true, male_rec_point, &m_rec_vector, &m_rec_pos);
				if(record_mf_rec_points==true){
					m_rec_points.push_back(male_rec_point);
					f_rec_points.push_back(female_rec_point);
				}
				if(record_m_rec_vector==true){
					write_rec_vector(m_rec_vector, m_rec_pos, m_rec_pos_array, m_rec_vector_out, o_genome_size);	
				}
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
	void write_rec_points(ofstream &m_rc_ofile, ofstream &f_rc_ofile, int rc_gen){
		//don't write generation 0, because it's confusing
		if(rc_gen!=0){
		//write female rec_points
			m_rc_ofile<<rc_gen<<"\t";
			for(int i_m_rp=0; i_m_rp<m_rec_points.size(); i_m_rp++){
				m_rc_ofile<<m_rec_points[i_m_rp]<<" ";
			}
			m_rc_ofile<<endl;
		//write male rec_points
			f_rc_ofile<<rc_gen<<"\t";
			for(int i_f_rp=0; i_f_rp<f_rec_points.size(); i_f_rp++){
				f_rc_ofile<<f_rec_points[i_f_rp]<<" ";
			}
			f_rc_ofile<<endl;
		}
	}
	void write_rec_vector(vector <double> *rv_rec_vector, vector <int> *rv_rec_points, vector <int> *rv_position_array, ofstream &rv_rec_vector_outfile, int rv_genome_size){

		if(rv_rec_points && rv_rec_vector){
			
			if(rv_position_array){
				int start_point=1;
				for(int i_rv_rec=0; i_rv_rec<rv_rec_vector->size()-1; i_rv_rec++){
					rv_rec_vector_outfile<<start_point<<"-"<<(*rv_position_array)[(*rv_rec_points)[i_rv_rec]]<<":"<<(*rv_rec_vector)[i_rv_rec]<<"\t";
					start_point=(*rv_position_array)[(*rv_rec_points)[i_rv_rec]]+1;
				}
				rv_rec_vector_outfile<<start_point<<"-"<<rv_genome_size<<":"<<(*rv_rec_vector)[rv_rec_vector->size()-1]<<endl;

			}
			else{
				for(int i_rv_rec=0; i_rv_rec<rv_rec_points->size(); i_rv_rec++){
					rv_rec_vector_outfile<<(*rv_rec_points)[i_rv_rec]<<":"<<(*rv_rec_vector)[i_rv_rec]<<"\t";
				}	
				rv_rec_vector_outfile<<endl;
			}
			
		}
	}
};

void read_config_file(string cfg_file_name, int &c_ngens, int &c_ninds, int &c_gen_size, double &c_mut_rate, double &c_mr_male_bias, int &c_sd_pos, int &c_rc_wind, double &c_min_ident, int &c_sstat_wind, int &c_sstat_gen, string &c_fst_ofile, string &c_dxy_ofile, string &c_dxx_ofile, double &c_rec_prob, double &c_var_effect, string &c_m_rp_ofile, string &c_f_rp_ofile, string &c_sdr_ofile, string &c_hs_conf_file, string &c_m_rv_ofile){
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
			else if(param=="rec_prob"){
				config_option>>c_rec_prob;
			}
			else if(param=="var_effect"){
				config_option>>c_var_effect;
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
			else if(param=="m_rp_ofile"){
				config_option>>c_m_rp_ofile;
			}
			else if(param=="f_rp_ofile"){
				config_option>>c_f_rp_ofile;
			}
			else if(param=="sdr_ofile"){
				config_option>>c_sdr_ofile;
			}
			else if(param=="hotspot_conf"){
				config_option>>c_hs_conf_file;
			}
			else if(param=="m_rv_ofile"){
				config_option>>c_m_rv_ofile;
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
	int wind_size=1000;
	//minimum identity to recombination prob of block to 0
	double min_ident=0.6;
	//genome size
	int gen_size=100000;
	//mutation rate
	double mutation_rate=2.5e-7;
	//male bias of mutation rate
	double mr_male_bias=1;
	//probability of crossing over landing in the simulated region
	double rec_prob=1;
	//relative effect of a single variant on recombination probabilities
	double var_effect=0;
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
	//male recombination positions file
	string m_rec_pos_ofile="m_rp.out";
	//female recombination positions file
	string f_rec_pos_ofile="f_rp.out";
	//sdr size ouput file
	string sdr_ofile="sdr_size.out";
	//name of the hotspot configuration file	
	string hotspot_conf_file="";

	string m_rv_ofile="";


	//run mode: 1=dynamic model(standard), 2=random recombination events, 3=hotspot model
	int run_mode=1;

	//read config file
	read_config_file(argv[1], ngens, n_inds, gen_size, mutation_rate, mr_male_bias, sd_loc, wind_size, min_ident, mean_rp_wind_size, sum_stat, fst_ofile, dxy_ofile, dxx_ofile, rec_prob, var_effect, m_rec_pos_ofile, f_rec_pos_ofile, sdr_ofile, hotspot_conf_file, m_rv_ofile);

	//open FST sumstat output
	ofstream fst_out;
	ofstream dxx_out;
	ofstream dxy_out;
	ofstream m_rp_out;
	ofstream f_rp_out;
	ofstream m_rv_out;
	//ofstream(sdr_out);
	dxx_out.open(dxx_ofile);
	dxy_out.open(dxy_ofile);
	fst_out.open(fst_ofile);
	m_rp_out.open(m_rec_pos_ofile);
	f_rp_out.open(f_rec_pos_ofile);
	if(m_rv_ofile!=""){
		m_rv_out.open(m_rv_ofile);
	}
	//sdr_out.open(sdr_ofile);
	fst_out<<"generation\tN_variants";
	dxx_out<<"generation";
	dxy_out<<"generation";
	m_rp_out<<"generation\trecombination positions"<<endl;
	f_rp_out<<"generation\trecombination positions"<<endl;
	//sdr_out<<"generation\tSDR_size"<<endl;
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
	if(hotspot_conf_file!=""){
		run_mode=3;
		initial_matrix.read_hotspot_conf(hotspot_conf_file);
	}
	population initial_pop(initial_sex_vector);
	for(int i_gen=0; i_gen<ngens; i_gen++){
		bool o_gen_rec_prob=false;
		bool o_gen_rec_m_rv=false;
		if(i_gen%sum_stat==0){
			o_gen_rec_prob=true;
			if(m_rv_ofile!=""){
				o_gen_rec_m_rv=true;
				m_rv_out<<"#gen "<<i_gen<<endl;
			}
		}

		initial_pop.make_offspring_gen(initial_matrix, wind_size, min_ident, r_generator, run_mode, n_inds, rec_prob, o_gen_rec_prob,m_rv_out,o_gen_rec_m_rv,var_effect, gen_size );
//		cout<<"Made offspring..."<<endl;
		initial_matrix.temp_remove_fixed_sites();
		initial_matrix.temp_add_mutations(mutation_rate, mr_male_bias, initial_pop.return_temp_female_pointer(), initial_pop.return_temp_male_pointer(), r_generator);
		initial_matrix.swap_arrays();
//		initial_matrix.print_array();
		if(i_gen%sum_stat==0){
			cout<<"Generation "<<i_gen<<endl;
			initial_matrix.write_window_fst(mean_rp_wind_size, initial_pop.return_male_pointer(), initial_pop.return_female_pointer(), i_gen, fst_out);
			initial_matrix.write_window_dxy(mean_rp_wind_size, initial_pop.return_female_pointer(), initial_pop.return_male_pointer(), dxx_out, dxy_out, i_gen);
			initial_pop.write_rec_points(m_rp_out,f_rp_out, i_gen);
			//sdr_out<<i_gen<<"\t";
			//initial_matrix.write_size_of_sdr(initial_pop.return_male_pointer(), sdr_out);
//			initial_matrix.print_position_array();
		}

	}
	fst_out.close();
	dxx_out.close();
	dxy_out.close();
	m_rp_out.close();
	f_rp_out.close();
	return 0;
};
