https://www.cheminventory.net



Admin Panel

User View



**Top Nav**

- Log In
- Search | Explore | Safety | My Chemicals | Order
- Search - structure, CAS number, name, synonyms, container, location, barcode, etc.
- Explore - list containers
- Safety - SDS
- My Chemicals - list chemicals
- Order - place orders, see recent orders, see current orders



**DB Design**

- Containers - see list of boxes, click on one and get list of chemicals
- My Chemicals - in my container, checked out by me, ordered by me
- Each Location (e.g. Xenon Pharma) has ID, name, address, building, floor, department_id(s)
- Each department has ID, name, function, contact information, shelf_id(s), user_id(s)
- Each shelf has ID, name, location description, box_id(s), deparment_id
- Each box has ID, name, bottle_id(s), shelf_id(s)
- Each bottle has ID, barcode, size, amount remaining, cost, supplier, order_id, box_id, chemical_id
- Each chemical has ID, name, CAS, structure, synonym(s), formula, description
- Each User has ID, role (head, chemist, biologist, safety, shipping/recieving, etc.), email, phone, department, container (optional), orders
- Each Order has ID, date placed, bottle_id(s)