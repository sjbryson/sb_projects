
"""
Samuel Joseph Bryson
Copyright 2026
"""

from neomodel import get_config, db



def test_connection():
    """
    Validates the connection to Neo4j and executes a simple Cypher query.
    """
    print("--- Neo4j Connection Test ---")
    config = get_config()
    config.database_url = 'bolt://neo4j:virus_db@localhost:7687'  # default
    # Check if the database_url is set
    if not config.database_url:
        print("Error: config.database_url is not set.")
        return
    print(f"Target URL: {config.database_url}")

    try:
        # run simple cypher query
        query = "RETURN 'Hello World' as message"
        results, meta = db.cypher_query(query)
        # neomodel returns results as a list of lists: [[ 'Hello World' ]]
        message = results[0][0]
        print(f"Success! Database responded with: {message}")
        print(f"Metadata columns: {meta}")
        
    except Exception as e:
        print(f"Connection Failed: {e}")

    finally:
        # Closes the entire driver/pool for the script
        db.close_connection()
        print(f"Connection Pool Closed.")

if __name__ == "__main__":
    #config = get_config()
    #config.database_url = 'bolt://neo4j:virus_db@localhost:7687'  # default
    test_connection()
    #db.close_connection()