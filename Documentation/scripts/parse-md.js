const path = require('path');
const fs = require('fs');
const docsDir = path.join(__dirname, '../docs');
const templatesDir = path.join(__dirname, '../src/templates');
const marked = require('marked');

for (const docPage of fs.readdirSync(docsDir)) {
  const docPath = path.join(docsDir, docPage);
  const doc = fs.readFileSync(docPath, 'utf-8');
  const parsed = marked(doc);
  const templateName = docPage.replace(/\..+/, '.js');
  const template = 'module.exports = `\n' + parsed + '\n`';
  // make template
  fs.writeFile(path.join(templatesDir, templateName), template, err => {
    if (err) {
      throw err;
    }
    console.log(`${templateName} created!`);
  });
}
