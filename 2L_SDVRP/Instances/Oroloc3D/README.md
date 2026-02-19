Test data for loading and routing 

Version: v1

The data contains 11 data sets. Each set is named following the format 'S{SEED_NUMBER}_{QUANTITY}T' where SEED_NUMBER is the random seed we used to sample the data and QUANTITY is the approximate total quantity in tonnes.

There are 3 sheets for each data set:
- `S{SEED_NUMBER}_{QUANTITY}T_items` - contains the input data.
    - `order_id` - integer id of order
    - `item_id` - integer id of item/package
    - `length`, `width`, `height` - size of item in millimeters, integer
    - `volume` - in m^3
    - `weight` - weight of item in kilograms
    - `delivery_address_id` - id of delivery address - the same for all packages in an order - see `addresses` sheet 
    
- `S{SEED_NUMBER}_{QUANTITY}T_sol_orders` - contains minimal info about solution at order level.
    - `order_id` - same as `order_id` in `S{SEED_NUMBER}_{QUANTITY}T_items`, used to identify order in solution
    - `truck_id` - id of truck in solution. Is -1 for undelivered orders.
    - `unload` - number of unload (1 is unloaded first at first destination, then 2 and so on). Is -1 for undelivered orders.
- `S{SEED_NUMBER}_{QUANTITY}T_sol_items` - contains detailed info about solution at item level.
    - `item_id`, `order_id`, - same as `item_id`, `order_id` in `S{SEED_NUMBER}_{QUANTITY}T_items`, used to identify item in solution (`item_id` should suffice)
    - `truck_id`, `unload` - same as `truck_id`, `unload` in `S{SEED_NUMBER}_{QUANTITY}T_sol_orders`
    - `x`, `y`, `z` - position of item inside trailer
    - `rotation_axis` - rotation of package around one axis when placed in trailer:
        - `0` - for rotation around x axis (swap width and height)
        - `1` - for rotation around y axis (swap length and height)
        - `2` - for rotation around z axis (swap length and width)
        - `-1` if not rotated

The sheet `addresses` contains a list of delivery addresses for all 11 data sets.
- `id` - id of address
- `country`, `city`, `country_code`, `delivery_address` - string info about address
- `latitude`, `longitude` coordinates of address

The first address with the coordinates (45.8454208059873, 24.9738267126969) was used as the start location.

The sheet `distances` contains road distances for all pairs of addresses.
- `id_1`, `id_2` - ids of addresses in the pair
- `distance` - road distance in meters


Parameters:
- Trailer size (length, width, height): (13500, 2450, 2700) millimeters
- Maximum loaded packages weight: 24000 kg
- Minimum loaded packages weight 21000 kg
- Max unloads: 4




