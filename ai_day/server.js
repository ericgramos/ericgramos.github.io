const express = require('express');
const fs = require('fs');
const path = require('path');
const app = express();

// Middleware to parse JSON bodies
app.use(express.json());

// Serve static files
app.use(express.static('.'));

// Handle registration submission
app.post('/submit-registration', (req, res) => {
    const { firstName, lastName, affiliation } = req.body;
    
    // Create CSV line
    const csvLine = `${firstName},${lastName},${affiliation}\n`;
    
    // Append to CSV file
    fs.appendFile('Participants.csv', csvLine, (err) => {
        if (err) {
            console.error('Error writing to CSV:', err);
            res.status(500).json({ error: 'Failed to save registration' });
            return;
        }
        res.status(200).json({ message: 'Registration successful' });
    });
});

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => {
    console.log(`Server running on port ${PORT}`);
}); 