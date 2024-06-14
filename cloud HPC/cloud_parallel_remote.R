## @knitr micro_cloud_parallel_remote_on_Linode_demo

# Introduction:----

# This script contains some examples for how a local R session can send code to
# and get it processed on a cloud HPC. This process is powered by the future 
# package. This functionality is showcased in one a future package documentation
# at https://furrr.futureverse.org/articles/remote-connections.html using AWS.
# Therefore, this script demos the same functionality on another cloud computing
# provider, Linode.

# Running code remotely on one Linode instances, no parallelisation:----

# 0. Example code: Fitting linear regression models:----

## clear R's global environment:
rm(list = ls())

## split the mtcars dataset by gear:
by_gear <- mtcars |>
  dplyr::group_split(gear) 
## fit the model using the 'purrr::map()' function (similar to lapply):
set.seed(1)
purrr_models <- purrr::map(
  .x = by_gear,
  .f = ~lm(mpg ~ cyl + hp + wt, data = .)
)
purrr_models                # print the models fitted using 'purrr::map()'
## fit the model using the 'furrr::future_map()' function:
future::plan(               # set the parallel processing plan as 'multisession'  
  strategy = future::multisession,
  workers = 2               # set the number of workers to parallel over
)
set.seed(1)
furrr_models <- furrr::future_map(
  .x = by_gear,
  .f = ~lm(mpg ~ cyl + hp + wt, data = .)
)
furrr_models                # print the models fitted using 'furrr::future_map()'

# 1. Get connection credentials:----
# 1.1. Save the Linode instance's IP address in a an object, "public_ip"

public_ip <- "176.58.117.117" # UIC_1_8GB

# 1.2. Save the path to your Linode's private SSH key, for example 'Linode_240514'

# check what ssh keys are set up in your system: 'ls -al ~/.ssh'
ssh_private_key_file <- "~/.ssh/Linode_240514"

# 2. Define a function that creates a connection to the remote cluster:----

connect_to_linode2 <- function(public_ip, ssh_private_key_file) {
  future::makeClusterPSOCK(
    
    # Public IP number of Linode instance
    workers = public_ip,
    
    # User name (could be 'ubuntu')
    user = "root",
    
    # Use private SSH key registered with Linode
    rshopts = c(
      "-o", "StrictHostKeyChecking=no",
      "-o", "IdentitiesOnly=yes",
      "-i", ssh_private_key_file
    ),
    
    rscript_args = c(
      # Set up .libPaths() for the 'ubuntu' user
      "-e", shQuote(paste0(
        "local({",
        "p <- Sys.getenv('R_LIBS_USER'); ",
        "dir.create(p, recursive = TRUE, showWarnings = FALSE); ",
        ".libPaths(p)",
        "})"
      )),
      # Install furrr
      "-e", shQuote("install.packages('furrr')")
    ),
    
    # Switch this to TRUE to see the code that is run on the workers without
    # making the connection
    dryrun = FALSE
  )
}
connect_to_linode <- function(public_ip, ssh_private_key_file) {
  future::makeClusterPSOCK(
    
    # Public IP number of Linode instance
    workers = public_ip,
    
    # User name (could be 'ubuntu')
    user = "root",
    
    # Use private SSH key registered with Linode
    rshopts = c(
      "-o", "StrictHostKeyChecking=no",
      "-o", "IdentitiesOnly=yes",
      "-i", ssh_private_key_file
    ),
    
    # Switch this to TRUE to see the code that is run on the workers without
    # making the connection
    dryrun = FALSE
  )
}

# 3. Create a connection to the cluster using the 'connect_to_linode()':---- 

cl <- connect_to_linode(public_ip, ssh_private_key_file)
# print the cluster object
cl

# 4. Run the code on the remote cluster:----
# 4.1. Define the parallel processing plan:

future::plan(                 
  strategy = future::cluster,   # set the parallel processing plan as 'cluster'
  workers = cl                  # pass the cluster connection defined above
)

# 4.2. Run the regression on the remote cluster:

furrr_models_cluster <- furrr::future_map(
  .x = by_gear, 
  .f = ~lm(mpg ~ cyl + hp + wt, data = .)
)
furrr_models_cluster

# 5. Revert back to a sequential plan and stop cluster:----

cl                                 # confirm cluster connection still active
parallel::stopCluster(cl = cl)     # stop cluster
cl                                 # confirm cluster connection stopped
future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential

# Running code remotely on one Linode instances, parallelisation:----

# 0. Run the test code sequentially:----
# 0.1. Ensure the plan is sequential:

future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential

# 0.2. Run the code locally in a sequentiol plan:

t1 <- proc.time()                  # keep track of execution starting time

res <- furrr::future_map(
  
  # Map over the two instances
  .x = c(1, 2), 
  
  .f = ~{
    
    outer_idx <- .x
    
    furrr::future_map_chr(
      
      # Each instance has 4 cores we can utilize
      .x = c(1, 2, 3, 4), 
      
      .f = ~{
        inner_idx <- .x
        Sys.sleep(2)
        paste0("Instance: ", outer_idx, " and core: ", inner_idx)
      }
    )
    
  }
)

t2 <- proc.time()                  # keep track of execution ending time
t2 - t1                            # print time it took to finish execution

# user  system elapsed 
# 0.073   0.119  16.061 

# 1. Create a connection to the cluster using the 'connect_to_linode()':---- 

cl <- connect_to_linode(public_ip, ssh_private_key_file)
# print the cluster object
cl

# 2. Run the code on the remote cluster:----
# 2.1. Define the parallel processing plan:

future::plan(
  strategy = list(
    # The outer plan instructs the outer future loop to distribute over the instances
    future::tweak(
      strategy = future::cluster,  # the outer future loop will be over clusters
      workers = cl                 # pass the cluster(s) connection defined above
    ), 
    
    # The inner plan instructs the inner future loop to run in parallel on each instance
    future::tweak(
      strategy = future::multisession, # the inner future loop will distribute over cores
      workers = 4                      # set the number of workers to parallel over
    )
  )
)

# 2.2. Run the code:

t1 <- proc.time()                  # keep track of execution starting time

res <- furrr::future_map(
  
  # Map over the two instances
  .x = c(1, 2), 
  
  .f = ~{
    
    outer_idx <- .x
    
    furrr::future_map_chr(
      
      # Each instance has 4 cores we can utilize
      .x = c(1, 2, 3, 4), 
      
      .f = ~{
        inner_idx <- .x
        Sys.sleep(2)
        paste0("Instance: ", outer_idx, " and core: ", inner_idx)
      }
    )
    
  }
)

t2 <- proc.time()                  # keep track of execution ending time
t2 - t1                            # print time it took to finish execution

#  user  system elapsed 
# 0.345   0.127   6.948

# 3. Revert back to a sequential plan and stop cluster:----

cl                                 # confirm cluster connection still active
parallel::stopCluster(cl = cl)     # stop cluster
cl                                 # confirm cluster connection stopped
future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential

# Running code remotely on multiple Linode instances, parallelisation:----

# 1. Create a connection to the cluster using the 'connect_to_linode()':---- 
# 1.1. add the ip of the second instance to the public_ip object:

public_ip                                    # print the current instance ip  
public_ip <- c(public_ip, "139.162.233.246") # UIC_1_16GB_8core
public_ip                                    # print the all instances ips

# 1.2. Create the remote connection to the instances:

# print the cluster object
cl
cl <- connect_to_linode(public_ip, ssh_private_key_file)
# print the cluster object
cl

# 2. Run the code on the remote cluster:----
# 2.1. Define the parallel processing plan:

future::plan(
  strategy = list(
    # The outer plan instructs the outer future loop to distribute over the instances
    future::tweak(
      strategy = future::cluster,  # the outer future loop will be over clusters
      workers = cl                 # pass the cluster(s) connection defined above
    ), 
    
    # The inner plan instructs the inner future loop to run in parallel on each instance
    future::tweak(
      strategy = future::multisession, # the inner future loop will distribute over cores
      workers = 4                      # set the number of workers to parallel over
    )
  )
)

# 2.2. Run the code:

t1 <- proc.time()                  # keep track of execution starting time

res <- furrr::future_map(
  
  # Map over the two instances
  .x = c(1, 2), 
  
  .f = ~{
    
    outer_idx <- .x
    
    furrr::future_map_chr(
      
      # Each instance has 4 cores we can utilize
      .x = c(1, 2, 3, 4), 
      
      .f = ~{
        inner_idx <- .x
        Sys.sleep(2)
        paste0("Instance: ", outer_idx, " and core: ", inner_idx)
      }
    )
    
  }
)

t2 <- proc.time()                  # keep track of execution ending time
t2 - t1                            # print time it took to finish execution

#  user  system elapsed 
# 0.245   0.102   4.550 

# 3. Revert back to a sequential plan and stop cluster:----

cl                                 # confirm cluster connection still active
parallel::stopCluster(cl = cl)     # stop cluster
cl                                 # confirm cluster connection stopped
future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential

# Running the HSSD PSA on Linode clusters:----

# 0. Source the Healthy-Sick-Sicker-Dead microsimulation model and PSA functions:----

source(
  file = file.path(here::here(), "extra_materials", "11_microsim", "cloud_parallel_remote_src.R")
)

# 1. Run and time the PSA locally, no clusters, but 4 workers:----
# 1.1. Set the parallel strategy to "multisession" with 4 workers:

future::plan(
  strategy = future::multisession, 
  workers = 4
)

# 1.2. Run and time the PSA:

system.time({
  df_psa_results <- run_psa_parallel(
    l_psa_parameters = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims = 16,
    num_workers = 4,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = FALSE,
    num_clusters = NULL
  )
})
#  user  system elapsed 
# 0.227   0.026   4.846 

# 2. Revert back to a sequential plan and stop cluster:----

future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential


# 3. Run the Healthy-Sick-Sicker-Dead microsimulation PSA in parallel over one linode clusters:----
# 3.1. Add one cluster to the 'public_ip' vector:

public_ip <- c(
  "139.162.233.246"  # UIC_1_16GB_8core
)

# 3.2. Create the cluster and connect to it:

# print the cluster object
if(exists("cl")) cl
cl <- connect_to_linode(public_ip, ssh_private_key_file)
# print the cluster object
cl

# 3.3. Define the parallel processing strategy:

future::plan(
  strategy = list(
    # The outer plan instructs the outer future loop to distribute over the instances
    future::tweak(
      strategy = future::cluster,  # the outer future loop will be over clusters
      workers = cl                 # pass the cluster(s) connection defined above
    ), 
    
    # The inner plan instructs the inner future loop to run in parallel on each instance
    future::tweak(
      strategy = future::multisession, # the inner future loop will distribute over cores
      workers = 4                      # set the number of workers to parallel over
    )
  )
)
future::plan()                    # print the current plan

# 3.4. Run PSA over the clusters and print the execution time:

system.time({
  df_psa_results_remote <- run_psa_parallel(
    l_psa_parameters = l_psa_parameters,
    n_sims = 16,
    num_workers = 4,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = TRUE,
    num_clusters = 1
  )
})
#  user  system elapsed 
# 0.753   0.068  15.754

# 3.5 Revert back to a sequential plan and stop cluster:

cl                                 # confirm cluster connection still active
parallel::stopCluster(cl = cl)     # stop cluster
cl                                 # confirm cluster connection stopped
future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential


# 4. Run the Healthy-Sick-Sicker-Dead microsimulation PSA in parallel over 4 linode clusters:----
# 4.1. Add all clusters to the 'public_ip' vector:

public_ip <- c(
  "176.58.117.117",  # UIC_1_8GB_4core
  "178.79.189.127",  # UIC_2_8GB_4core
  "178.79.169.76",   # UIC_3_8GB_4core
  "139.162.233.246"  # UIC_1_16GB_8core
)

# 4.2. Create the cluster and connect to it:

# print the cluster object
if(exists("cl")) cl
cl <- connect_to_linode(public_ip, ssh_private_key_file)
# print the cluster object
cl

# 4.3. Define the parallel processing strategy - multisession:

future::plan(
  strategy = list(
    # The outer plan instructs the outer future loop to distribute over the instances
    future::tweak(
      strategy = future::cluster,  # the outer future loop will be over clusters
      workers = cl                 # pass the cluster(s) connection defined above
    ), 
    
    # The inner plan instructs the inner future loop to run in parallel on each instance
    future::tweak(
      strategy = future::multisession, # the inner future loop will distribute over cores
      workers = 4                      # set the number of workers to parallel over
    )
  )
)
future::plan()                    # print the current plan

# 4.4. Run PSA over the clusters (multisession) and print the execution time:

system.time({
  df_psa_results_remote2 <- run_psa_parallel(
    l_psa_parameters = l_psa_parameters,
    n_sims = 16,
    num_workers = 4,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = TRUE,
    num_clusters = 4
  )
})
#  user  system elapsed 
# 0.788   0.095  13.608 

# 4.5. Define the parallel processing strategy - multicore:

future::plan(
  strategy = list(
    # The outer plan instructs the outer future loop to distribute over the instances
    future::tweak(
      strategy = future::cluster,  # the outer future loop will be over clusters
      workers = cl                 # pass the cluster(s) connection defined above
    ), 
    
    # The inner plan instructs the inner future loop to run in parallel on each instance
    future::tweak(
      strategy = future::multicore, # the inner future loop will distribute over cores
      workers = 4                   # set the number of workers to parallel over
    )
  )
)
future::plan()                    # print the current plan

# 4.6. Run PSA over the clusters (multicore) and print the execution time:

system.time({
  df_psa_results_remote3 <- run_psa_parallel(
    l_psa_parameters = l_psa_parameters,
    n_sims = 16,
    num_workers = 4,
    parallel_method = "multicore",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = TRUE,
    num_clusters = 4
  )
})
#  user  system elapsed 
# 0.634   0.094   9.662  

# 4.7. Compare local and remote PSA outputs:

testthat::expect_equal(df_psa_results, df_psa_results_remote)
testthat::expect_equal(df_psa_results, df_psa_results_remote2)
testthat::expect_equal(df_psa_results, df_psa_results_remote3)

# 4.8. Run 100K individuals in 10K PSA on clusters and print the execution time:

system.time({
  df_psa_results_remote4 <- run_psa_parallel(
    l_psa_parameters = l_psa_parameters,
    n_sims = 1e4,
    num_workers = 4,
    parallel_method = "multicore",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = TRUE,
    num_clusters = 4
  )
})
# 1e2
#  user  system elapsed 
# 2.382   0.194  35.994
# 1e3
#   user  system elapsed 
# 20.361   1.431 301.169 
  
# 4.5. Revert back to a sequential plan and stop cluster

cl                                 # confirm cluster connection still active
parallel::stopCluster(cl = cl)     # stop cluster
cl                                 # confirm cluster connection stopped
future::plan(
  strategy = future::sequential
)
future::plan()                     # confirm the plan is back to sequential
