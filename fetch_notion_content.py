import os
from notion_client import Client

# Initialize the Notion client with the API key
notion = Client(auth=os.getenv('NOTION_API_KEY'))

# Replace this with your actual page ID
page_id = os.getenv('NOTION_PAGE_ID')

# Fetch the page content
page = notion.pages.retrieve(page_id=page_id)

# Fetch the blocks of the page (content)
blocks = notion.blocks.children.list(block_id=page_id)

# Parse content into Markdown
markdown_content = []

for block in blocks['results']:
    if block['type'] == 'paragraph':
        markdown_content.append(block['paragraph']['text'][0]['plain_text'])
    elif block['type'] == 'heading_1':
        markdown_content.append(f"# {block['heading_1']['text'][0]['plain_text']}")
    elif block['type'] == 'heading_2':
        markdown_content.append(f"## {block['heading_2']['text'][0]['plain_text']}")
    elif block['type'] == 'heading_3':
        markdown_content.append(f"### {block['heading_3']['text'][0]['plain_text']}")
    elif block['type'] == 'bulleted_list_item':
        markdown_content.append(f"- {block['bulleted_list_item']['text'][0]['plain_text']}")
    elif block['type'] == 'numbered_list_item':
        markdown_content.append(f"1. {block['numbered_list_item']['text'][0]['plain_text']}")
    # Add more block types as needed

# Write content to a Markdown file
with open('notion-sync/notion_page.md', 'w') as file:
    file.write('\n'.join(markdown_content))
