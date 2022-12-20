PGG_data_setup <- function(path, V23, V126, V128, V138){
  raw <- read.csv(path, 
                  skip = 3, header = T)
  
  ### removing punishment condition 
  raw <- raw %>% 
    filter(p == "N-experiment")
  
  ### all cities in experiments - and corresponding country
  cities = unique(raw$city)
  nations = c("Russia", "Belarus", "USA", "Oman", "Switzerland", "Denmark", 
              "England", "Ukraine", "Saudi Arabia", "Tyrkey", "China", "South Korea",
              "Germany", "Greece", "Australia", "Switzerland")
  
  ### cities we could find nationa data for
  city = c("Samara", "Boston", "St. Gallen", "Zurich", "Nottingham", "Dnipropetrovs'k", 
           "Istanbul", "Chengdu", "Seoul", "Bonn", "Melbourne")
  # nation ID
  nation = c(1 , 2 , 3 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10)
  id = seq(1,10,1)
  
  ### merging data 
  nation_data <- data.frame(id, V23, V126, V128, V138) %>% 
    rename(nation = id)
  city_df = data.frame(city, nation)
  nation_data = merge(nation_data, city_df)
  
  # making final df - automatically removes countries without data
  df <- merge(raw, nation_data , on = "city") %>% 
    # select every third row
    filter(row_number() %% 3 == 1)
  
  return (df)
}

calculate_winnings_df <- function(df, pi, Nation_df){
  # sum of contributions for each group on each trial times pi
  win_gs <- df %>% 
    group_by(period, groupid) %>% 
    summarise(sum_con_trial = sum(senderscontribution)) %>% 
    mutate(win_trial = sum_con_trial*pi)
  
  # sum of winnings for each group
  winnings_m <- win_gs %>% 
    group_by(groupid) %>% 
    summarise(win = sum(win_trial))
  
  # merging with nation to later find winnings at nation-level
  winnings <- merge(winnings_m, Nation_df)
  
  return(winnings)
  
}
